import os
import urllib.request

from damply import dirs

# Downloads the file and saves it locally


def main(
	ccl_meta_url: str = 'https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/cellinfo_beta.txt',
	ccl_file: str = 'cellinfo_beta.txt',
	cpd_meta_url: str = 'https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/compoundinfo_beta.txt',
	cpd_file: str = 'compoundinfo_beta.txt',
	gene_meta_url: str = 'https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/geneinfo_beta.txt',
	gene_file: str = 'geneinfo_beta.txt',
	sig_meta_url: str = 'https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/siginfo_beta.txt',
	sig_file: str = 'siginfo_beta.txt',
	# inst_meta_url:str = "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/instinfo_beta.txt",
	# inst_file:str = "instinfo_beta.txt",
	cpd_sig_url: str = 'https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/level5/level5_beta_ctl_n58022x12328.gctx',
	data_file: str = 'level5_beta_ctl_n58022x12328.gctx',
):

	os.system('Rscript fetch_adb.R')
	for name, url, destination in zip(
		[ccl_meta_url, cpd_meta_url, gene_meta_url, sig_meta_url, cpd_sig_url],
		[ccl_file, cpd_file, gene_file, sig_file, data_file],
	):
		urllib.request.urlretrieve(url, dirs.RAWDATA / destination)


if __name__ == '__main__':
	main()
