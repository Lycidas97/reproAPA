"""Join the sealed raw-source checks with final contributor checks once available."""
import json

from verify_APA_delivery import HERE, SOURCE
from verify_marker_delivery import LOCK, sha
from verify_unique_support import verify_unique_support


raw = json.loads((HERE / "raw_delivery_validation.json").read_text())
assert raw["status"] == "PASS"
assert raw["verifier_sha256"] == sha(HERE / "verify_APA_delivery.py")
assert raw["inputs_sha256"] == sha(SOURCE / "inputs.npz")
classification = json.loads((HERE / "delivery_validation.json").read_text())
assert classification["classification_check"]["status"] == "PASS"
acceptance = json.loads((SOURCE / "analysis_validation.json").read_text())
assert acceptance["status"] == "COMPLETE_AWAITING_INDEPENDENT_QA"
assert raw["raw_profile_checks"]["source_profile_rows"] == acceptance["n_group_profile_rows"]
assert raw["complete_gene_eligibility"]["eligibility_rows"] == acceptance["n_gene_pair_qualification_rows"]
assert sha(SOURCE / "PROTOCOL.zh.md") == LOCK
assert sha(SOURCE / "ADDENDUM.zh.md") == "ddbfe768ef6210ebda71cdf584e70dd36a0f5e9f3691ed48d5716445758b956a"
weights = json.loads((HERE / "frozen_weights_validation.json").read_text())
assert weights["status"] == "PASS"
for name, digest in weights["source_sha256"].items():
    assert sha(SOURCE / name) == digest
assert weights["measurable_pairs"] == acceptance["n_measurable_pairs"]
assert weights["gradient_pairs"] == acceptance["n_gradient_pairs"]
vectors = json.loads((HERE / "frozen_pair_vectors_validation.json").read_text())
assert vectors["status"] == "PASS" and vectors["pair_rows"] == acceptance["n_measurable_pairs"]
gradient_reference = json.loads((HERE / "gradient_reference_validation.json").read_text())
assert gradient_reference["status"] == "PASS" and gradient_reference["gradient_pairs"] == acceptance["n_gradient_pairs"]
overlap = json.loads((HERE / "shared_MC_joint_support_validation.json").read_text())
assert overlap["n_eligible_pairs_lost_after_distinct_MC_dedup"] == 15
unique_support = verify_unique_support()
report = dict(raw, status="PASS", contributor_weights=weights,
              full_pair_vector_checks=vectors,
              gradient_fixed_original_reference_check=gradient_reference,
              marker_classification=classification["classification_check"],
              full_primary_axes=classification["full_primary_axes"],
              frozen_qualification_rule="group-contribution numeric replay; distinct-source-MC support claim reserved for unique_MC_support sensitivity",
              shared_source_MC_diagnostic=overlap, unique_source_MC_qualification=unique_support,
              source_analyses=acceptance,
              source_sha256={name: sha(SOURCE / name) for name in [
                  "PROTOCOL.zh.md", "ADDENDUM.zh.md", "cohort_coverage.tsv", "gene_background.tsv.gz",
                  "primary_axis_features.tsv.gz", "library_apa_eligibility.tsv", "pair_eligibility.tsv.gz",
                  "gene_group_profiles.tsv.gz", "pair_effects.tsv.gz", "gradient_pair_effects.tsv.gz",
                  "bridge_pair_weights.tsv.gz", "analysis_validation.json", "QA-ADDENDUM.zh.md",
                  "unique_support_pair_weights.tsv.gz", "unique_support_route_effects.tsv.gz",
                  "unique_support_bridge_route_effects.tsv.gz"]},
              finalizer_sha256=sha(__file__))
(HERE / "apa_delivery_validation.json").write_text(json.dumps(report, indent=2) + "\n")
classification.update(phase="independent_marker_raw_APA_weights_PASS", production_APA_checked=True,
                      independent_APA_report="apa_delivery_validation.json", independent_scope_complete=True,
                      complete_delivery_verified=True,
                      verification_scope="all classification/cohort/axis/pre-gene/matrix keys and contributor weights; fixed stratified raw PAS/profile/pair qualification reconstruction; route summaries reviewed separately by root")
(HERE / "delivery_validation.json").write_text(json.dumps(classification, indent=2) + "\n")
print(json.dumps(report, indent=2))
