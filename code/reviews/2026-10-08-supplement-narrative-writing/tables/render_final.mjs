import fs from 'node:fs/promises';
import path from 'node:path';
import {FileBlob, SpreadsheetFile} from '@oai/artifact-tool';
process.on('uncaughtException',e=>{console.error(e?.message??String(e));process.exitCode=1;});
process.on('unhandledRejection',e=>{console.error(e?.message??String(e));process.exitCode=1;});
const here=path.dirname(new URL(import.meta.url).pathname);
const wb=await SpreadsheetFile.importXlsx(await FileBlob.load(path.join(here,'final-skim.xlsx')));
const ranges={S9a:'A1:J19',S9b:'A1:J19',S9c:'A1:J19',S9d:'A1:J19',S9e:'A1:J19',Contents:'A56:E62',Sources:'A159:F163',p50:'A21:L34',decomposition:'I8:V18'};
const inspections=[];
for(const [name,range] of Object.entries(ranges)){
 const sheetName=['p50','decomposition'].includes(name)?'S9e':name;
 const img=await wb.render({sheetName,range,scale:1.4,format:'png'});
 await fs.writeFile(path.join(here,`${name}-final-preview.png`),new Uint8Array(await img.arrayBuffer()));
 const check=await wb.inspect({kind:'table',range:`${sheetName}!${range}`,include:'values,formulas',tableMaxRows:20,tableMaxCols:15,maxChars:2200});
 inspections.push({name,sheetName,range,ndjson:check.ndjson});
}
await fs.writeFile(path.join(here,'final-artifact-inspect.json'),JSON.stringify(inspections,null,2)+'\n');
console.log('Rendered final workbook read-only samples.');
