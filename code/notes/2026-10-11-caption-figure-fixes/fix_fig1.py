"""Figure 1: 1d legend '1306 orthologs' -> '1,306 orthologs'; 1c colourbar '0.03' -> '+0.03' (house sign style)."""
import json
from fixlib import Job
J = Job(1)
J.replace("1306", (300, 572, 365, 584), "1,306 orthologs", note="1d legend")
J.drop_word("orthologs", (300, 572, 365, 584), note="redrawn with the count")
J.replace('0.03', (185, 569, 210, 581), '+0.03', align='center', note='1c colourbar end label')
rec = J.run()
print(json.dumps(rec, indent=1, ensure_ascii=False))
json.dump(rec, open('out/receipt.fig1.json', 'w'), indent=1, ensure_ascii=False)
