# E57 local runbook

Required local parent results are E54 and E56. Place the preregistered E57 script in the repository scripts directory.

First run the script with the self-test flag. The expected marker is E57_SYNTHETIC_SIX_REPLICA_SELF_TEST_PASS.

Then run the same script with the run flag and capture stdout/stderr to e57_ninemock_six48k.log.

Only R5/R6 are newly evaluated. R1-R4 are reused from E56. The result/checkpoint path is source_data/e57_conditioned_full_eligible_ninemock_six48k_result.json. The script checkpoints after every expensive term and can be resumed with the same run command.
