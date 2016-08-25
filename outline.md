# Microbial Innovations

### Overview:

--

#### A. Geological Signals of Innovation

1. Representative Trees

	- [New Tree of Life](http://www.nature.com/articles/nmicrobiol201648)
	- [Silva Living Tree Project](https://www.arb-silva.de/projects/living-tree/)
	- [RDP Release 11.4](https://rdp.cme.msu.edu/misc/rel10info.jsp)
	- Tree must be ultrametric (treePL?)

	- 22-25 Aug 2016: Fasta sequences for trees were downloaded.
	- Silva LTP includes non-uniuqe record IDs; I wrote a python script (check_unique.py) to parse the records in a multi-FASTA file and determine if record IDs are unique. If non-unique record IDs exist, the script saves only the first record.

	- FastTree being used to generate phylogenetic trees

			+ parameters: -gtr -nt -gamma

			`FastTree -gtr -nt -gamma LTPs123_unique.fasta > LTPs123.tree `
			`Optimize all lengths: LogLk = -2028364.001 Time 1747.74
			Gamma(20) LogLk = -2030537.456 alpha = 0.361 rescaling lengths by 7.156   
			Total time: 2119.25 seconds Unique: 11901/11933 Bad splits: 81/11898 Worst delta-LogLk 12.421`
			`FastTree -gtr -nt -gamma nmicrobiol201648-s7.txt > NTL.tree`
			`Optimize all lengths: LogLk = -812454.853 Time 111.90
			Gamma(20) LogLk = -816606.980 alpha = 0.487 rescaling lengths by 1.161s   
			 Total time: 126.69 seconds Unique: 1854/1871 Bad splits: 0/1851`

	 - treePL being used to make trees ultrametric

2. Tree Patterns

	- Branch Lengths
	- Bursts: size & distribution

3. Geological Time Events

4. Do "burts" map onto/near reasonable geological time events?

5. Do present-day traits agree with Burst-Geology inferences?

--

#### B. Microbial Innovations

1. Create the new, recalibrated tree

2. Are "bursts" a proxy for functional groups?

--

### Outcomes:
