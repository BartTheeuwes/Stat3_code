for f in *.fq.gz; do i="${f%.fq.gz}"; mv "$f" "${i//./_}.fastq.gz"; done


rename _s_1_r_2. _S1_L001_R2_001. *.gz
rename _s_1_r_1. _S1_L001_R1_001. *.gz
rename _s_1_i_1. _S1_L001_I2_001. *.gz
rename _s_1_i_2. _S1_L001_I2_001. *.gz

rename _s_2_r_2. _S2_L001_R2_001. *.gz
rename _s_2_r_1. _S2_L001_R1_001. *.gz
rename _s_2_i_1. _S2_L001_I2_001. *.gz
rename _s_2_i_2. _S2_L001_I2_001. *.gz