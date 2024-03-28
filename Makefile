setup:
	@rm -rf bin/
	@bash scripts/download_third_party.sh
	@conda env update
	@bash scripts/post_configure_conda_env.sh 

download_test_data:
	@bash scripts/download_test_data.sh

run_test:
	@bash scripts/run_test.sh

download_paper_data:
	@bash scripts/download_paper_data.sh

run_paper_analysis:
	@bash scripts/run_paper_analysis.sh