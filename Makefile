setup:
	@rm -rf bin/
	@bash scripts/download_third_party.sh
	@conda env update
	@bash scripts/post_configure_conda_env.sh