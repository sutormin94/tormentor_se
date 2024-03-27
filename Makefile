setup:
	@rm -rf bin/
	@bash scripts/download_third_party.sh
	@conda env update
	@bash scripts/post_configure_conda_env.sh

download_test_data:
	@cd tests & fasterq-dump  

test:
	@tormentor --reads tests/reads_1.fastq tests/reads_2.fastq -o tests/results/