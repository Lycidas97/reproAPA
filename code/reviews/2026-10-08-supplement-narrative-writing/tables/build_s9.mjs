import fs from 'node:fs/promises';
import path from 'node:path';
import {FileBlob, SpreadsheetFile, Workbook} from '@oai/artifact-tool';
process.on('uncaughtException',e=>{console.error(e?.message??String(e));process.exitCode=1;});
process.on('unhandledRejection',e=>{console.error(e?.message??String(e));process.exitCode=1;});
const here=path.dirname(new URL(import.meta.url).pathname);
const payload=JSON.parse(await fs.readFile(path.join(here,'s9-payload.json'),'utf8'));
const reference=await SpreadsheetFile.importXlsx(await FileBlob.load(path.join(here,'reference-skim.xlsx')));
for(const name of ['Contents','S4e']) {
  const blob=await reference.render({sheetName:name,range:'A1:H17',scale:1.4,format:'png'});
  await fs.writeFile(path.join(here,`reference-${name}.png`),new Uint8Array(await blob.arrayBuffer()));
}
const wb=Workbook.create();
const font='Arial', ink='#1F2A36', accent='#23415E', muted='#6B7785', band='#F2F5F8';
const log=[];
function col(n){let s='';for(n++;n;n=Math.floor((n-1)/26))s=String.fromCharCode(65+(n-1)%26)+s;return s;}
function typed(v){
  if(v==='')return null;
  if(typeof v==='number')return Number.isFinite(v)?v:String(v);
  if(typeof v==='boolean')return String(v);
  if(/^[-+]?\d+$/.test(v) && (BigInt(v)>BigInt(Number.MAX_SAFE_INTEGER)||BigInt(v)<BigInt(Number.MIN_SAFE_INTEGER)))return v;
  if(/^[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?$/.test(v)) {
    const x=Number(v);if(Number.isFinite(x))return x;
  }
  return v;
}
for(const [name,title] of Object.entries(payload.titles)){
  const sh=wb.worksheets.add(name);sh.showGridLines=false;sh.tabColor=accent;
  sh.getRange('A1:H3').format.font={name:font,size:9,color:muted};
  sh.getRange('A1').values=[[`Table ${name}. ${title}`]];
  sh.getRange('A1').format.font={name:font,size:13,bold:true,color:accent};
  sh.getRange('A2').values=[['Table S9. Analyses specified after the initial results. '+payload.analyses[name]+'.']];
  sh.getRange('A3').values=[['Initial results preceded these analyses. Missing values remain NA; no statistical results were recalculated.']];
  sh.freezePanes.freezeRows(3);sh.freezePanes.freezeColumns(1);
  let r=5,maxCols=8;
  const bs=payload.blocks.filter(b=>b.sheet===name);
  for(const b of bs){
    const width=b.columns.length,end=col(width-1);maxCols=Math.max(maxCols,width);
    sh.getRange(`A${r}`).values=[[`${b.block}  ${b.title} (${b.rows.toLocaleString('en-US')} rows)`]];
    sh.getRange(`A${r}`).format.font={name:font,size:10,bold:true,color:accent};r++;
    if(b.notes){sh.getRange(`A${r}`).values=[[b.notes]];sh.getRange(`A${r}`).format.font={name:font,size:9,italic:true,color:muted};r++;}
    if(b.filter!=='All rows; no statistical recomputation.'){
      // Exact field projections are documented in Sources; short row conditions stay beside the data.
      const short=b.filter.split(' Columns: ')[0];
      sh.getRange(`A${r}`).values=[['Rows shown: '+short]];
      sh.getRange(`A${r}`).format.font={name:font,size:9,italic:true,color:muted};r++;
    }
    const header=r,range=`A${header}:${end}${header+b.rows}`;
    const values=[b.columns,...b.values.map(row=>row.map(typed))];
    sh.getRange(range).values=values;
    sh.getRange(range).format.font={name:font,size:10,color:ink};
    sh.getRange(range).format.verticalAlignment='top';
    sh.getRange(`A${header}:${end}${header}`).format={fill:accent,font:{name:font,size:10,bold:true,color:'#FFFFFF'},wrapText:true,verticalAlignment:'center',rowHeight:42};
    for(let j=0;j<width;j++){
      const letter=col(j),body=sh.getRange(`${letter}${header+1}:${letter}${header+b.rows}`);
      const vals=values.slice(1).map(row=>row[j]).filter(v=>v!==null&&v!=='NA');
      const numeric=vals.length>0 && vals.every(v=>typeof v==='number');
      body.format.horizontalAlignment=numeric?'right':'left';
      if(numeric){
        const ints=vals.every(v=>Number.isInteger(v)&&Math.abs(v)<1e15);
        body.setNumberFormat(ints?'0':'0.0000');
        // Preserve readable small p values without changing stored precision.
        if(/(?:^|_)(?:p|q|raw_p|p_value|q_BH|probability)(?:_|$)/i.test(b.columns[j])||b.columns[j].includes('probability'))body.setNumberFormat('0.000E+00');
      }
      const textWidth=Math.min(40,Math.max(9,b.columns[j].length*.82+3));
      sh.getRange(`${letter}1:${letter}1`).format.columnWidth=textWidth;
    }
    sh.getRange(`A${header+1}:${end}${header+b.rows}`).conditionalFormats.addCustom(`MOD(ROW()-${header},2)=0`,{fill:band});
    const tbl=sh.tables.add(range,true,`T_${name}_${b.block.split('.')[1]}`);
    tbl.showFilterButton=true;tbl.showBandedRows=false;
    log.push({sheet:name,block:b.block,header_row:header,rows:b.rows,columns:width,source:b.source,table_name:tbl.name});
    r=header+b.rows+3;
  }
  // Wide notes and titles span empty cells; source tables retain separate fields.
  sh.getRange(`A1:${col(maxCols-1)}3`).format.rowHeight=20;
  console.log(`Authored ${name}: ${bs.length} blocks`);
}
wb.recalculate();
const inspected=[];
for(const name of Object.keys(payload.titles)){
  const check=await wb.inspect({kind:'table',range:`${name}!A1:J18`,include:'values,formulas',tableMaxRows:18,tableMaxCols:10,maxChars:1800});
  inspected.push({sheet:name,ndjson:check.ndjson});
  const img=await wb.render({sheetName:name,range:'A1:J19',scale:1.4,format:'png'});
  await fs.writeFile(path.join(here,`${name}-preview.png`),new Uint8Array(await img.arrayBuffer()));
}
const errors=await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:100},maxChars:1500});
const out=await SpreadsheetFile.exportXlsx(wb);await out.save(path.join(here,'S9-native.xlsx'));
await fs.writeFile(path.join(here,'s9-blocks.json'),JSON.stringify(log,null,2)+'\n');
await fs.writeFile(path.join(here,'artifact-inspect.json'),JSON.stringify({inspected,formulaErrors:errors.ndjson},null,2)+'\n');
console.log(JSON.stringify({sheets:Object.keys(payload.titles),blocks:log.length,output:path.join(here,'S9-native.xlsx')}));
