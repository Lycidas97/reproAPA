import fs from 'node:fs/promises';
import path from 'node:path';
import assert from 'node:assert/strict';
import { Workbook, SpreadsheetFile } from '@oai/artifact-tool';

// The statistical adapters supply typed values; no statistics are recomputed in Excel.
const root = path.resolve(process.argv[2]);
const payload = JSON.parse(await fs.readFile(path.join(root, 'workbook-inputs.json'), 'utf8'));
const wb = Workbook.create();
const previews = path.join(root, 'workbook-previews');
await fs.mkdir(previews, { recursive: true });
const audit = [];
const ink = '#17324D';
const historical = '#985D20';
function literal(value) {
  if (value !== null && typeof value === 'object') return JSON.stringify(value);
  if (typeof value === 'string' && value.startsWith('=')) return `'${value}`;
  return value;
}
function addTable(name, columns, data, scope) {
  assert(columns.length > 0 && data.every(row => row.length === columns.length));
  const s = wb.worksheets.add(name);
  s.showGridLines = false;
  s.freezePanes.freezeRows(1);
  s.tabColor = scope.includes('historical') ? historical : ink;
  s.getRangeByIndexes(0, 0, 1, columns.length).values = [columns];
  for (let start = 0; start < data.length; start += 2000) {
    const block = data.slice(start, start + 2000).map(row => row.map(literal));
    s.getRangeByIndexes(start + 1, 0, block.length, columns.length).values = block;
  }
  const used = s.getRangeByIndexes(0, 0, data.length + 1, columns.length);
  used.format.font = { name: 'Arial', size: 10, color: '#202D38' };
  used.format.rowHeight = 17;
  used.format.columnWidth = 18;
  const widths = columns.map(key => /source|path|origin|file|gene_symbol|name|scope|status|reason/i.test(key) ? 32 : 18);
  const headerHeight = Math.max(42, ...columns.map((key, j) => 14 * Math.ceil(key.length / widths[j]) + 8));
  s.getRangeByIndexes(0, 0, 1, columns.length).format = {
    fill: scope.includes('historical') ? historical : ink,
    font: { name: 'Arial', size: 10, bold: true, color: '#FFFFFF' },
    wrapText: true, rowHeight: headerHeight, verticalAlignment: 'center', horizontalAlignment: 'center',
  };
  for (let j = 0; j < columns.length; j++) {
    const key = columns[j];
    const width = widths[j];
    s.getRangeByIndexes(0, j, data.length + 1, 1).format.columnWidth = width;
    const sample = data.map(row => row[j]).filter(value => value !== null);
    if (sample.length && sample.every(value => typeof value === 'number')) {
      const tiny = sample.some(value => value !== 0 && Math.abs(value) < 0.0001);
      const format = sample.every(Number.isInteger) ? '#,##0' : tiny || /(^p_|^q_|_p$|_q$|pval|p_value|q_value|empirical_p|q_BH)/i.test(key) ? '0.0000E+00' : '0.0000';
      if (data.length) s.getRangeByIndexes(1, j, data.length, 1).setNumberFormat(format);
    }
  }
  audit.push({ name, columns: columns.length, rows: data.length, scope });
  return s;
}
const notes = [
  ['Run', payload.run_id],
  ['Input fingerprint', payload.input_fingerprint],
  ['S1–S7', 'Frozen upstream GEX metadata; unchanged.'],
  ['S8–S30', 'Recomputed PAS analyses. The Index lists every source part and its evidence scope.'],
  ['Historical audit', payload.historical_scope],
  ['Missing values', 'Blank numeric cells mean unavailable or inapplicable. They are not zeros or nonsignificant tests. Read each status column.'],
  ['Precision', payload.full_precision],
  ['Large tables', 'Parts marked complete_TSV_only are delivered losslessly as TSV. S29-statistical-source-data.zip contains primary numerical records; S29-execution-provenance.zip contains run configurations and analysis code. S18-v3-reference-control-source-data.zip contains the separately recomputed reference-control statistics.'],
  ['Statistics', 'P values, adjusted P values, null ranges and effect intervals are calculated by the registered analysis code. XLSX contains exported results.'],
];
const readme = addTable('Read me', ['Field', 'Definition'], notes, 'documentation');
readme.getRange('A1:A10').format.columnWidth = 24;
readme.getRange('B1:B10').format.columnWidth = 100;
readme.getRange('B2:B10').format.wrapText = true;
readme.getRange('A2:B10').format.rowHeight = 44;
const index = addTable('Index', ['Table', 'Part', 'Delivery', 'Evidence scope', 'Source', 'Data file', 'Note'],
  payload.sheets.map(s => [`S${s.table}`, s.name, s.role, s.evidence_scope, s.origin, s.data_file, s.reason || '']), 'source_index');
if (payload.sheets.length) {
  index.getRangeByIndexes(1, 0, payload.sheets.length, 7).format.wrapText = true;
  index.getRangeByIndexes(1, 0, payload.sheets.length, 7).format.rowHeight = 56;
}
for (const entry of payload.sheets) {
  if (entry.role === 'complete_TSV_only') continue;
  assert.equal(entry.role, 'complete_table');
  const value = JSON.parse(await fs.readFile(path.join(root, entry.source), 'utf8'));
  addTable(entry.name, value.columns, value.data, entry.evidence_scope);
}
wb.recalculate();
const errors = await wb.inspect({ kind: 'match', searchTerm: '#REF!|#DIV/0!|#VALUE!|#NAME\\?|#NUM!|#NULL!|#SPILL!|#CALC!',
  options: { useRegex: true, maxResults: 300 }, maxChars: 12000 });
await fs.writeFile(path.join(root, 'workbook-error-scan.json'), JSON.stringify(errors, null, 2));
// Render a bounded view of every worksheet; full tables remain present in the workbook.
for (const item of audit) {
  const inspection = await wb.inspect({ kind: 'region', sheetId: item.name, range: 'A1:H10', maxChars: 3500, tableMaxRows: 10, tableMaxCols: 8 });
  await fs.writeFile(path.join(previews, `${item.name}.inspection.json`), JSON.stringify(inspection, null, 2));
  const preview = await wb.render({ sheetName: item.name, range: 'A1:H10', scale: 1.4, format: 'png' });
  await fs.writeFile(path.join(previews, `${item.name}.png`), new Uint8Array(await preview.arrayBuffer()));
}
const destination = path.join(root, 'table', 'reproAPA-polyApipe-S1-S30.xlsx');
await fs.mkdir(path.dirname(destination), { recursive: true });
await (await SpreadsheetFile.exportXlsx(wb)).save(destination);
await fs.writeFile(path.join(root, 'workbook-export.json'), JSON.stringify({ run_id: payload.run_id,
  input_fingerprint: payload.input_fingerprint, output: destination, sheets: audit,
  rendered_every_sheet: true, visual_review_complete: false }, null, 2));
console.log(JSON.stringify({ output: destination, worksheets: audit.length, visual_review_required: true }));
