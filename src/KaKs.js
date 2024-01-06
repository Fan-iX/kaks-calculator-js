/************************************************************
* Filename: KaKs.js
* Abstract: Declaration of KAKS class including several methods.

* Version: js-2.0
* Author: Fanix
* Date: December.11, 2022

* Adapted from: KaKs_Calculator/KaKs.cpp and KaKs_Calculator/KaKs.h
* Copyright (C) 2009, Beijing Institute of Genomics of CAS
* All rights reserved.

* Modified Version: 2.0.1
* Modified Author: Kristian K Ullrich
* Modified Date: April.29, 2020

* Version: 2.0
* Author: Da-Peng Wang(wangdp@big.ac.cn), Yu-Bin Zhang (ybzhang@big.ac.cn)
* Date: Jun.1, 2009

* Version: 1.0
* Author: Zhang Zhang (zhanghzhang@genomics.org.cn)
* Date: Jan.21, 2005

*************************************************************/

class KAKS extends Base {
	constructor() {
		super();
		this.items = ["Sequence", "Method", "Ka", "Ks", "Ka/Ks",
			"P-Value(Fisher)", "Length", "S-Sites", "N-Sites", "Fold-Sites(0:2:4)",
			"Substitutions", "S-Substitutions", "N-Substitutions", "Fold-S-Substitutions(0:2:4)", "Fold-N-Substitutions(0:2:4)",
			"Divergence-Time", "Substitution-Rate-Ratio(rTC:rAG:rTA:rCG:rTG:rCA/rCA)", "GC(1:2:3)", "ML-Score", "AICc",
			"Akaike-Weight", "Model"];
		this.method_name = [
			"NG", "LWL", "LPB", "MLWL", "MLPB", "GY", "YN", "MYN", "MS", "MA", "GNG", "GLWL", "GLPB", "GMLWL", "GMLPB", "GYN", "GMYN"
		];
		this.method_ref = [
			"Nei, M. and Gojobori, T. (1986) Mol. Biol. Evol., 3, 418-426.",
			"Li, W.H., Wu, C.I. and Luo, C.C. (1985) Mol. Biol. Evol., 2, 150-174.",
			"Li, W.H. (1993) J. Mol. Evol., 36, 96-99.    Pamilo, P. and Bianchi, N.O. (1993) Mol. Biol. Evol., 10, 271-281.",
			"Tzeng, Y.H., Pan, R. and Li, W.H. (2004) Mol. Biol. Evol., 21, 2290-2298.",
			"Tzeng, Y.H., Pan, R. and Li, W.H. (2004) Mol. Biol. Evol., 21, 2290-2298.",
			"Goldman, N. and Yang, Z. (1994) Mol. Biol. Evol., 11, 725-736.",
			"Yang, Z. and Nielsen, R. (2000) Mol. Biol. Evol., 17, 32-43.",
			"Zhang, Z., Li, J. and Yu, J. (2006) BMC Evolutionary Biology, 6, 44.",
			"Model Selection according to the AICc",
			"Model Averaging on a set of candidate models",
			"Wang, DP., Zhang, S., He, FH., Zhu, J.,Hu, SN. and Yu, J.(2009) Genomics, Proteomics and Bioinformatics. In press.",
			"Wang, DP., Zhang, S., He, FH., Zhu, J.,Hu, SN. and Yu, J.(2009) Genomics, Proteomics and Bioinformatics. In press.",
			"Wang, DP., Zhang, S., He, FH., Zhu, J.,Hu, SN. and Yu, J.(2009) Genomics, Proteomics and Bioinformatics. In press.",
			"Wang, DP., Zhang, S., He, FH., Zhu, J.,Hu, SN. and Yu, J.(2009) Genomics, Proteomics and Bioinformatics. In press.",
			"Wang, DP., Zhang, S., He, FH., Zhu, J.,Hu, SN. and Yu, J.(2009) Genomics, Proteomics and Bioinformatics. In press.",
			"Wang, DP., Zhang, S., He, FH., Zhu, J.,Hu, SN. and Yu, J.(2009) Genomics, Proteomics and Bioinformatics. In press.",
			"Wang, DP., Wan, HL., Zhang, S. and Yu, J. (2009) Biology Direct, 4:20 (16 June 2009)"
		];
		this.tempt = 0;
		this.Initialize();
	}

	Initialize() {
		this.none = this.ng86 = this.gng86 = this.lpb93 = this.glpb93 = this.lwl85 = this.glwl85 = this.mlwl85 = this.gmlwl85 = this.mlpb93 = this.gmlpb93 = this.yn00 = this.gyn00 = this.gy94 = this.myn = this.gmyn = this.ms = this.ma = false;
		this.number = 0;
		this.details = "";
		this.result = this.seq1 = this.seq2 = "";
		seq_name = ""
		genetic_code = 1;
	}

	/* Get GCC of entire sequences GC[0] and of three codon positions GC[1,2,3] */
	getGCContent(str) {
		let i, j;
		GC.fill(0);
		for (i = 0; i < str.length; i += 3) {
			let codon = str.substr(i, 3);
			for (j = 0; j < 3; j++) {
				if (codon[j] == 'G' || codon[j] == 'C') GC[j + 1]++;
			}
		}
		GC[0] = this.sumArray(GC, 4, 1) / str.length;
		for (i = 1; i < 4; i++) GC[i] /= (str.length / 3);
	}

	/****************************************************
	* Function: ReadCalculateSeq
	* Input Parameter: string
	* Output: Read sequences, check sequences' validity
			  and calculate Ka and Ks.
	* Return Value: True if succeed, otherwise false.
	 
	* Note: Using axt file for low memory
	*****************************************************/
	ReadCalculateSeq(data) {
		this.showParaInfo();
		this.result = this.getTitleInfo();
		for (let seqs of data.trim().split(/\n{2,}/)) {
			let name = seqs.substr(0, seqs.indexOf('\n'));
			let str = seqs.substr(seqs.indexOf('\n') + 1).replace(/\n/g, "");
			this.getGCContent(str);
			this.checkValid(name, str.substr(0, str.length / 2), str.substr(str.length / 2, str.length / 2));
			this.calculateKaKs();
			console.log(name + " [OK]");
		}
	}

	/**************************************************
	* Function: checkValid
	* Input Parameter: string, string, string
	* Output: Check validity of pairwise sequences
	* Return Value: True if succeed, otherwise false. 
	***************************************************/
	checkValid(name, str1, str2) {
		let i;

		//Check whether (sequence length)/3==0
		if (str1.length != str1.length || str1.length % 3 != 0 || str2.length % 3 != 0) {
			throw "Error. The size of two sequences in '" + name + "' is not equal.";
		}

		//Delete gap and stop codon
		let arr1 = str1.toUpperCase().match(/.{3}/g)
		let arr2 = str2.toUpperCase().match(/.{3}/g)
		for (i = 0; i < arr1.length; i++) {
			if (
				/[^ATCG]/g.test(arr1[i]) |
				/[^ATCG]/g.test(arr2[i]) |
				this.getAminoAcid(arr1[i]) == "!" |
				this.getAminoAcid(arr2[i]) == "!"
			) {
				arr1[i] = arr2[i] = ""
			}
		}

		//pass value into private variables
		this.seq1 = arr1.join("");
		this.seq2 = arr2.join("");
		//pass value into extern variables
		seq_name = name;
		length = this.seq1.length;
	}

	/**************************************************
	* Function: Run
	* Input Parameter: string, object
	* Output: Calculate Ka/Ks and output.
	* Return Value: void 
	***************************************************/
	Run(data, argv) {

		let time_start = new Date()
		this.parseParameter(argv)
		this.ReadCalculateSeq(data)
		//Time used for running
		let t = (new Date() - time_start) / 1000 | 0;
		let h = t / 3600 | 0, m = (t % 3600) / 60 | 0, s = t - (t / 60 | 0) * 60;

		//Print on display
		console.log("Mission accomplished. (Time elapsed:" + (h ? h + ":" : "") + m + ":" + s + ")");
		return this.result
	}


	/**************************************************
	* Function: parseParameter
	* Input Parameter: int, const char* []
	* Output: Parse the input parameters
	* Return Value: bool 
	***************************************************/
	parseParameter(argv) {
		//No parameter
		if (argv.length == 0) throw 1;
		else if (argv.length == 1) {//help information
			let temp = argv[0];
			if (temp == "-H" || temp == "-h") {
				return helpInfo();
			}
			else {
				throw 1;
			}
		}
		else {
			//parse parameters
			let codeflag = 0;
			for (let i = 0; i < argv.length; i++) {

				let temp = argv[i].toUpperCase();
				//Genetic Code Table
				if (temp == "-C") {
					if ((i + 1) < argv.length && codeflag == 0) {
						genetic_code = +argv[++i];
						if (genetic_code < 1 || genetic_code > NCODE || transl_table[2 * (genetic_code - 1)].length < 1)
							throw 1;
						codeflag++;
					}
					else {
						throw 1;
					}
				}
				//Algorithm(s) selected
				else if (temp == "-M") {
					if ((i + 1) > argv.length) throw 1;
					temp = argv[++i].toUpperCase();

					if (temp == "NONE") this.none = true;
					else if (temp == "NG") this.ng86 = true;
					else if (temp == "GNG") this.gng86 = true; //added by zhangyubin
					else if (temp == "LWL") this.lwl85 = true;
					else if (temp == "GLWL") this.glwl85 = true;
					else if (temp == "LPB") this.lpb93 = true;
					else if (temp == "GLPB") this.glpb93 = true;
					else if (temp == "MLPB") this.mlpb93 = true;
					else if (temp == "GMLPB") this.gmlpb93 = true;
					else if (temp == "MLWL") this.mlwl85 = true;
					else if (temp == "GMLWL") this.gmlwl85 = true;
					else if (temp == "GY") this.gy94 = true;
					else if (temp == "YN") this.yn00 = true;
					else if (temp == "GYN") this.gyn00 = true;
					else if (temp == "MYN") this.myn = true;
					else if (temp == "GMYN") this.gmyn = true;
					else if (temp == "MS") this.ms = true;
					else if (temp == "MA") this.ma = true;
					else if (temp == "ALL") {
						this.ng86 = this.gng86 = this.lpb93 = this.glpb93 = this.lwl85 = this.glwl85 = this.mlwl85 = this.gmlwl85 = this.mlpb93 = this.gmlpb93 = this.gy94 = this.yn00 = this.gyn00 = this.myn = this.gmyn = this.ms = this.ma = true; //zhangyubin  added
					}
					else throw 1;
				}
				else throw 1;
			}

			//Default: use ma to to calculate Ka and Ks
			if (!(this.none + this.ng86 + this.gng86 + this.lpb93 + this.glpb93 + this.lwl85 + this.glwl85 + this.mlwl85 + this.gmlwl85 + this.mlpb93 + this.gmlpb93 + this.gy94 + this.yn00 + this.gyn00 + this.myn + this.gmyn + this.ms + this.ma)) { //zhangyubin added
				this.ma = true;
			}
		}
	}
	calculateKaKs() {
		//Estimate Ka and Ks
		if (this.none) this.start_NONE(0);
		if (this.ng86 || this.gng86) this.start_NG86(0);
		if (this.lwl85 || this.glwl85) this.start_LWL85(0);
		if (this.mlwl85 || this.gmlwl85) this.start_MLWL85(0);
		if (this.lpb93 || this.glpb93) this.start_LPB93(0);
		if (this.mlpb93 || this.gmlpb93) this.start_MLPB93(0);
		if (this.gy94) this.start_GY94(0);
		if (this.yn00 || this.gyn00) this.start_YN00(0);
		if (this.myn || this.gmyn) this.start_MYN(0);
		if (this.ms || this.ma) this.start_MSMA(0);

		if (this.gmyn || this.gng86 || this.gyn00 || this.glwl85 || this.gmlwl85 || this.glpb93 || this.gmlpb93) {

			let fka, fks, fkaks = 0, temp = 0;
			let tresult;
			tresult = this.result;
			if (temp == 0) {
				tresult = tresult.substr(tresult.indexOf('\n') + 1);
			}
			temp++;

			//	while ((tresult.indexOf('\n'))>0 && tresult.length()>0 && tresult!="") {	
			//	for (int num=0;num<6;num++)
			//	{

			let j, k;
			j = tresult.indexOf('\n');
			let linecontent = tresult.substr(0, j + 1);
			linecontent += '\t';
			tresult = tresult.substr(j + 1);

			k = 1;
			while ((j = linecontent.indexOf('\t')) > 0) {
				let temp = linecontent.substr(0, j);

				if (k == 3) {
					fka = +temp;
				}
				if (k == 4) {
					fks = +temp;
				}
				if (k == 5) {
					fkaks = +temp;
					//goto att;
				}

				linecontent = linecontent.substr(j + 1);
				k++;
			}
			//	}

			if ((!this.myn) && (!this.ng86) && (!this.lwl85) && (!this.mlwl85) && (!this.lpb93) && (!this.mlpb93) && (!this.yn00) && (!this.gy94) && (!this.ms) && (!this.ma)) {
				if (this.tempt == 0) {
					// this.result = this.getTitleInfo();//for Gamma choices result delete the first line 
					//Linux version don't need this line
					this.tempt++;
				} else {
					// this.result = "";//for Gamma choices result delete the first line 
				}
				if (this.gy94) this.start_GY94(-1);
				if (this.ms || this.ma) this.start_MSMA(-1);

			}

			if (fkaks < 1) {
				//	result="";
				//Estimate Ka and Ks
				if (this.none) this.start_NONE(-1);
				if (this.gng86) this.start_NG86(-1);
				if (this.glwl85) this.start_LWL85(-1);
				if (this.gmlwl85) this.start_MLWL85(4);
				if (this.glpb93) this.start_LPB93(1);
				if (this.gmlpb93) this.start_MLPB93(1);
				//	if (this.gy94)	this.start_GY94(0);
				if (this.gyn00) this.start_YN00(4);
				if (this.gmyn) this.start_MYN(20);

				//	if (this.ms||this.ma)	this.start_MSMA(0);
			} else if (fkaks > 1) {

				//	result="";
				//Estimate Ka and Ks
				if (this.none) this.start_NONE(-1);
				if (this.gng86) this.start_NG86(6);
				if (this.glwl85) this.start_LWL85(0.2);
				if (this.gmlwl85) this.start_MLWL85(0.6);
				if (this.glpb93) this.start_LPB93(1);
				if (this.gmlpb93) this.start_MLPB93(1);
				//	if (this.gy94)	this.start_GY94(0);
				if (this.gyn00) this.start_YN00(-1);
				if (this.gmyn) this.start_MYN(-1);

				//	if (this.ms||this.ma)	this.start_MSMA(0);
			} else if (fkaks == 1) {
				if (this.none) this.start_NONE(-1);
				if (this.gng86) this.start_NG86(-1);
				if (this.glwl85) this.start_LWL85(-1);
				if (this.gmlwl85) this.start_MLWL85(-1);
				if (this.glpb93) this.start_LPB93(-1);
				if (this.gmlpb93) this.start_MLPB93(-1);
				//	if (this.gy94)	this.start_GY94(-1);
				if (this.gyn00) this.start_YN00(-1);
				if (this.gmyn) this.start_MYN(-1);
				//	if (this.ms||this.ma)	this.start_MSMA(-1);
			}

		}
	}
	start_NONE(GAMMA) {
		let zz = new NONE()
		zz.GAMMA = GAMMA
		this.result += zz.Run(this.seq1, this.seq2)
	}
	start_NG86(GAMMA) {
		let zz = new NG86()
		zz.GAMMA = GAMMA
		this.result += zz.Run(this.seq1, this.seq2)
	}
	start_LWL85(GAMMA) {
		let zz = new LWL85()
		zz.GAMMA = GAMMA
		this.result += zz.Run(this.seq1, this.seq2)
	}
	start_MLWL85(GAMMA) {
		let zz = new MLWL85()
		zz.GAMMA = GAMMA
		this.result += zz.Run(this.seq1, this.seq2)
	}
	start_LPB93(GAMMA) {
		let zz = new LPB93()
		zz.GAMMA = GAMMA
		this.result += zz.Run(this.seq1, this.seq2)
	}
	start_MLPB93(GAMMA) {
		let zz = new MLPB93()
		zz.GAMMA = GAMMA
		this.result += zz.Run(this.seq1, this.seq2)
	}
	start_GY94(GAMMA) {
		let zz = new GY94("HKY")
		zz.GAMMA = GAMMA
		this.result += zz.Run(this.seq1, this.seq2)
	}
	start_YN00(GAMMA) {
		let zz = new YN00()
		zz.GAMMA = GAMMA
		this.result += zz.Run(this.seq1, this.seq2)
	}
	start_MYN(GAMMA) {
		let zz = new MYN()
		zz.GAMMA = GAMMA
		this.result += zz.Run(this.seq1, this.seq2)
	}
	/************************************************
	* Function: start_MSMA
	* Input Parameter: void
	* Output: Calculate Ka and Ks using the method of 
			  model selection or model averaging.
	* Return Value: void
	*************************************************/
	start_MSMA(GAMMA) {

		let result4MA = [];	//generated by MS and used by MA

		//Model Selection
		let zz1 = new MS();
		zz1.GAMMA = GAMMA;
		let tmp = zz1.Run(this.seq1, this.seq2, result4MA, this.details);
		if (this.ms) {
			this.result += tmp;
		}

		//Model Averaging
		if (this.ma) {
			let zz2 = new MA();
			zz2.GAMMA = GAMMA;
			this.result += zz2.Run(this.seq1, this.seq2, result4MA);
		}
	}

	/************************************************
	* Function: showParaInfo
	* Input Parameter: 
	* Output: print on display
	* Return Value: 
	*************************************************/
	showParaInfo() {
		let info = [];
		let methods = ["Method(s):"];
		if (this.none) methods.push("NONE");
		if (this.ng86) methods.push("NG");
		if (this.gng86) methods.push("GNG");
		if (this.lwl85) methods.push("LWL");
		if (this.glwl85) methods.push("GLWL");
		if (this.mlwl85) methods.push("MLWL");
		if (this.gmlwl85) methods.push("GMLWL");

		if (this.lpb93) methods.push("LPB");
		if (this.glpb93) methods.push("GLPB");
		if (this.mlpb93) methods.push("MLPB");
		if (this.gmlpb93) methods.push("GMLPB");
		if (this.gy94) methods.push("GY");
		if (this.yn00) methods.push("YN");
		if (this.gyn00) methods.push("GYN");
		if (this.myn) methods.push("MYN");
		if (this.gmyn) methods.push("GMYN");	//zhangyubin added
		if (this.ms) methods.push("MS");
		if (this.ma) methods.push("MA");
		info.push(methods.join(" "));
		info.push("Genetic code:" + transl_table[2 * (genetic_code - 1) + 1])
		info.push("Please wait while reading sequences and calculating...")
	}

	/************************************************
	* Function: getTitleInfo
	* Input Parameter: 
	* Output: get title information of outputing file
	* Return Value: string
	*************************************************/
	getTitleInfo() {
		return this.items.join("\t") + "\n";
	}

}
