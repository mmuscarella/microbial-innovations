# Microbial Innovations

### Overview:

--

#### A. Geological Signals of Innovation

1. Representative Trees

	- [New Tree of Life](http://www.nature.com/articles/nmicrobiol201648)
	- [Silva Living Tree Project](https://www.arb-silva.de/projects/living-tree/)
	- [RDP Release 11.4 Type Strains](https://rdp.cme.msu.edu/misc/rel10info.jsp)
	- Tree must be ultrametric (treePL?)

	- 22-25 Aug 2016: Fasta sequences for trees were downloaded.
	- RDP Notes: RDP reference FASTA dataset with Type stains from cultured and uncultured sources with size > 1200 bp and good quality (25 Aug 2016)
	- Silva LTP includes non-uniuqe record IDs; I wrote a python script (check_unique.py) to parse the records in a multi-FASTA file and determine if record IDs are unique. If non-unique record IDs exist, the script saves only the first record.

	- FastTree being used to generate phylogenetic trees

			+ parameters: -gtr -nt -gamma

			`FastTree -gtr -nt -gamma LTPs123_unique.fasta > LTPs123.tree`
				`ML-NNI round 10: LogLk = -2028381.534 NNIs 157 max delta 9.60 Time 1615.08 (final)
				Optimize all lengths: LogLk = -2028364.161 Time 1688.31
				Gamma(20) LogLk = -2030536.709 alpha = 0.361 rescaling lengths by 7.156   
				Total time: 2058.36 seconds Unique: 11901/11933 Bad splits: 82/11898 Worst delta-LogLk 12.425`

			`FastTree -gtr -nt -gamma nmicrobiol201648-s7.txt > NTL.tree`
				`ML-NNI round 9: LogLk = -812498.450 NNIs 18 max delta 30.35 Time 109.14 (final)
				Optimize all lengths: LogLk = -812494.339 Time 112.14
				Gamma(20) LogLk = -816625.111 alpha = 0.487 rescaling lengths by 1.161s   
				Total time: 125.47 seconds Unique: 1854/1871 Bad splits: 5/1851 Worst delta-LogLk 7.870`

			`FastTree -gtr -nt -gamma rdp_download_9752seqs.fa > RDP.tree`
				`ML-NNI round 9: LogLk = -1754062.900 NNIs 137 max delta 12.05 Time 1029.89 (final)
				Optimize all lengths: LogLk = -1754049.823 Time 1063.57
				Gamma(20) LogLk = -1757170.448 alpha = 0.295 rescaling lengths by 4.795   
				Total time: 1261.53 seconds Unique: 9501/9752 Bad splits: 79/9498 Worst delta-LogLk 8.486`

			`classify.seqs(fasta=LTPs123_unique.fasta, template=trainset14_032015.pds.fasta, taxonomy=trainset14_032015.pds.tax, cutoff=80, probs=T)

			classify.seqs(fasta=LTPs123_unique.fasta, template=silva.nr_v123.align, taxonomy=silva.nr_v123.tax, cutoff=50, probs=F)

			classify.seqs(fasta=nmicrobiol201648_s7.txt, template=trainset14_032015.pds.fasta, taxonomy=trainset14_032015.pds.tax)


			classify.seqs(fasta=rdp_download_9752seqs.fa, template=trainset14_032015.pds.fasta, taxonomy=trainset14_032015.pds.tax)`

	 - treePL being used to make trees ultrametric


	 		+ treePL requires the number of sites. I'm going to assume this would be the number of bases. So, I need a python script that determines the number of bases in each record. This should be the same because these are alignments, but we will see.

			-

		+ We have had to make changes to the treePL code. The original code ignores branches that are under a set threshold. However, we do not want to ignore the distance represented by these brances. So we rewrote the code to add the distance of branches which are removed to the downstream branches at the given node.

		+ We have now had some discussion about identifying bursts in ultrametric trees.  

2. Tree Patterns

	- Branch Lengths

		+ branch length distribution (see how James did this)

	- Bursts: size & distribution

		+

3. Geological Time Events

4. Do "burts" map onto/near reasonable geological time events?

5. Do present-day traits agree with Burst-Geology inferences?

--

#### B. Microbial Innovations

1. Create the new, recalibrated tree

2. Are "bursts" a proxy for functional groups?

--

### Outcomes:
