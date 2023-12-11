/*********************************************************
* Filename: LWL85.js
* Abstract: Definition of LWL85 and Modified LWL85 class.

* Version: js-2.0
* Author: Fanix
* Date: December.11, 2022

* Adapted from: KaKs_Calculator/LWL85.cpp and KaKs_Calculator/LWL85.h
* Copyright (C) 2009, Beijing Institute of Genomics of CAS
* All rights reserved.

* Modified Version: 2.0.1
* Modified Author: Kristian K Ullrich
* Modified Date: April.29, 2020

* Version: 2.0
* Author: Da-Peng Wang(wangdp@big.ac.cn), Yu-Bin Zhang (ybzhang@big.ac.cn)
* Date: Jun.1, 2009

* Version: 1.0
* Author: Zhang Zhang (zhang.zhang@yale.edu)
* Date: Feb., 2005

  References: 
  Li WH, Wu CI, Luo CC  (1985)  A new method for
  estimating synonymous and nonsynonymous rates of nucleotide 
  substitution considering the relative likelihood of nucleotide
  and codon changes. Mol. Biol. Evol. 2:150-174.

  Tzeng Y-H, Pan R, Li W-H  (2004)  Comparison of Three Methods
  for Estimating Rates of Synonymous and Nonsynonymous Nucleotide
  Substitutions. Mol. Biol. Evol. 21:2290-2298.
**********************************************************/

class LWL85 extends Base {
	constructor() {
		super()
		this.name = "LWL"
		this.K.fill(0)
		this.A = new Array(5).fill(0)
		this.B = new Array(5).fill(0)
		this.Pi = new Array(5).fill(0)
		this.Qi = new Array(5).fill(0)
		this.Si_temp = new Array(5).fill(0)
		this.Vi_temp = new Array(5).fill(0)
	}

	/************************************************
	* Function: CountSiteAndDiff
	* Input Parameter: codon1, codon2
	* Output: Calculate synonymous and nonsynonymous
			 sites and differences between two codons.
	* Return Value: void
	*************************************************/
	CountSiteAndDiff(codon1, codon2) {

		let i, j, k;
		let num = 0, diff = new Array(5).fill(0);
		let temp1 = "", temp2 = "";
		let sum = 1, stop = 0;

		//Count sites
		for (i = 0; i < CODONLENGTH; i++) {
			this.L[this.getCodonClass(codon1, i)]++;
			this.L[this.getCodonClass(codon2, i)]++;
		}

		//Count differences
		for (i = 0; i < 3; i++) {
			diff[i] = -1;
			if (codon1[i] != codon2[i]) {
				diff[num] = i;
				num++;
			}
		}

		//Two codons are identical
		if (num == 0) return;

		//Sum is the number of evolution pathway.
		for (i = 1; i <= num; i++) sum *= i;

		this.snp += num;

		for (i = 0; i < CODONLENGTH; i++)
			this.Si_temp[2 * i] = this.Vi_temp[2 * i] = 0.0;

		//One difference
		if (num == 1) {
			this.TransitionTransversion(codon1, codon2, diff[0]);
		}

		//Two differences: 2 pathways for evolution (i,j)
		if (num == 2) {
			for (i = 0; i < num; i++) {
				for (j = 0; j < num; j++) {
					if (i != j) {
						temp1 = codon1.split("");
						temp1[diff[i]] = codon2[diff[i]];
						temp1 = temp1.join("")
						if (this.getAminoAcid(temp1) != '!') {
							//pathway: codon1 <-> temp1 <-> codon2
							this.TransitionTransversion(codon1, temp1, diff[i]);
							this.TransitionTransversion(temp1, codon2, diff[j]);
						}
						else {
							stop++;
						}
					}
				}
			}
		}

		//Three differences: 6 pathways for evolution (i,j,k)
		if (num == 3) {
			for (i = 0; i < 3; i++) {
				for (j = 0; j < 3; j++) {
					for (k = 0; k < 3; k++) {
						if ((i != j) && (i != k) && (j != k)) {
							temp1 = codon1.split("");
							temp1[diff[i]] = codon2[diff[i]];
							temp2 = temp1.map(x => x);
							temp2[diff[j]] = codon2[diff[j]];
							temp1 = temp1.join("")
							temp2 = temp2.join("")
							if (this.getAminoAcid(temp1) != '!' && this.getAminoAcid(temp2) != '!') {
								//pathway: codon1 <-> temp1 <-> temp2 <-> codon2
								this.TransitionTransversion(codon1, temp1, diff[i]);
								this.TransitionTransversion(temp1, temp2, diff[j]);
								this.TransitionTransversion(temp2, codon2, diff[k]);
							}
							else
								stop++;
						}
					}
				}
			}
		}

		//Add pair-codon's differences to Si and Vi
		for (i = 0; i < CODONLENGTH && (sum - stop) > 0; i++) {
			this.Si[2 * i] += (this.Si_temp[2 * i] / (sum - stop));
			this.Vi[2 * i] += (this.Vi_temp[2 * i] / (sum - stop));
		}
	}

	/************************************************
	* Function: getCodonClass
	* Input Parameter: codon, position(0,1,2)
	* Output: return 0,2,4-fold of codon at a given position.
	* Return Value: 0 or 2 or 4
	*************************************************/
	getCodonClass(codon, pos) {
		let i;
		let codonClass = 0;

		for (i = 0; i < 4; i++) {
			if (i != this.convertChar(codon[pos])) {
				let temp = codon.split("");
				temp[pos] = this.convertInt(i);
				temp = temp.join("")
				if (this.getAminoAcid(temp) != '!' && this.getAminoAcid(temp) == this.getAminoAcid(codon)) {
					codonClass++;
				}
			}
		}
		if (codonClass > 0 && codonClass < 3) {
			codonClass = 2;
		}
		else if (codonClass == 3) {
			codonClass = 4;
		}

		return codonClass;
	}

	/************************************************
	* Function: TransitionTransversion
	* Input Parameter: codon1, codon2, position(0,1,2)
	* Output: Calculate synonymous and nonsynonymous differences
			of two codons at a given position.
	* Return Value: int
	* Note: Follow kakstools.py sent from Prof.Li
	*************************************************/
	TransitionTransversion(codon1, codon2, pos) {

		//1:synonymous, 0:nonsynonymous
		let isSyn = 0;

		/* Follow kakstools.py sent from Prof.Li WH */
		//CGA, CGG, AGA, AGG
		if ((codon1 == "CGA" && codon2 == "AGA" && pos == 0) || (codon1 == "AGA" && codon2 == "CGA" && pos == 0)) {
			isSyn = 1;
		}
		if ((codon1 == "CGG" && codon2 == "AGG" && pos == 0) || (codon1 == "AGG" && codon2 == "CGG" && pos == 0)) {
			isSyn = 1;
		}
		if (isSyn == 1) {
			let c1 = this.getCodonClass(codon1, pos);
			let c2 = this.getCodonClass(codon2, pos);
			this.Si_temp[c1] += 0.5;	//different from the following 
			this.Vi_temp[c2] += 0.5;

			return 1;
		}

		/* Normal situation */
		//Synonymous: T(0)<->C(1), A(2)<->G(3) 	
		let sum = this.convertChar(codon1[pos]) + this.convertChar(codon2[pos]);
		if (sum == 5 || sum == 1)
			isSyn = 1;

		let class1 = this.getCodonClass(codon1, pos);
		let class2 = this.getCodonClass(codon2, pos);
		if (isSyn == 1) {
			this.Si_temp[class1] += 0.5;
			this.Si_temp[class2] += 0.5;
		}
		if (isSyn == 0) {
			this.Vi_temp[class1] += 0.5;
			this.Vi_temp[class2] += 0.5;
		}

		return 1;
	}

	/************************************************
	* Function: preProcess
	* Input Parameter: seq1, seq2
	* Output: preprocess for Run
	* Return Value: void
	*************************************************/
	preProcess(seq1, seq2) {

		let i;
		let ts = 0, tv = 0;
		let ai = new Array(5).fill(0), bi = new Array(5).fill(0);

		for (i = 0; i < seq1.length; i += 3) {
			this.CountSiteAndDiff(seq1.substr(i, 3), seq2.substr(i, 3));
		}

		for (i = 0; i < 5; i += 2) {

			ai[i] = bi[i] = 0.0;

			ts += this.Si[i];
			tv += this.Vi[i];

			this.L[i] = this.L[i] / 2.0;
			this.Pi[i] = this.Si[i] / this.L[i];
			this.Qi[i] = this.Vi[i] / this.L[i];

			if ((1 - 2 * this.Pi[i] - this.Qi[i]) > 0 && (1 - 2 * this.Qi[i]) > 0) {
				ai[i] = 1 / (1 - 2 * this.Pi[i] - this.Qi[i]);


				if ((Math.abs(this.GAMMA - 0.2) < SMALLVALUE) || this.GAMMA == -1) {
					this.name = "GLWL";

				}
				if ((Math.abs(this.GAMMA - 0.6) < SMALLVALUE) || this.GAMMA == 4 || this.GAMMA == -2) {
					this.name = "GMLWL";

				}
				if ((this.GAMMA == 1 || this.GAMMA == -3) && (this.name == "LPB")) {
					this.name = "GLPB";

				}
				if ((this.GAMMA == 1 || this.GAMMA == -4) && (this.name == "MLPB")) {
					this.name = "GMLPB";

				}

				if ((Math.abs(this.GAMMA - 0.2) < SMALLVALUE) || (Math.abs(this.GAMMA - 0.6) < SMALLVALUE) || this.GAMMA == 4 || this.GAMMA == 1) {

					ai[i] = 1 - 2 * this.Pi[i] - this.Qi[i];
					ai[i] = Math.pow(ai[i], -1.0 / this.GAMMA) - 1;
				}
				bi[i] = 1 / (1 - 2 * this.Qi[i]);

				if ((Math.abs(this.GAMMA - 0.2) < SMALLVALUE) || (Math.abs(this.GAMMA - 0.6) < SMALLVALUE) || this.GAMMA == 4 || this.GAMMA == 1) {
					bi[i] = 1 - 2 * this.Qi[i];
					bi[i] = Math.pow(bi[i], -1.0 / this.GAMMA) - 1;
				}
				//zhangyubin add 
			}

			if ((Math.abs(this.GAMMA - 0.2) < SMALLVALUE) || (Math.abs(this.GAMMA - 0.6) < SMALLVALUE) || this.GAMMA == 4 || this.GAMMA == 1) {
				this.B[i] = this.GAMMA * 0.5 * bi[i];
				this.A[i] = this.GAMMA * 0.5 * ai[i] - 0.25 * this.GAMMA * bi[i];
				this.K[i] = this.A[i] + this.B[i];
			}
			else {
				if (ai[i] > 0 && bi[i] > 0) {
					if (Math.log(bi[i]) >= 0) {
						this.B[i] = 0.5 * Math.log(bi[i]);
					}

					if ((0.5 * Math.log(ai[i]) - 0.25 * Math.log(bi[i])) >= 0) {
						this.A[i] = 0.5 * Math.log(ai[i]) - 0.25 * Math.log(bi[i]);
					}
					this.K[i] = this.A[i] + this.B[i];
				}
			}
		}

		if (tv > 0)
			this.kappa = 2 * ts / tv;
		else
			this.kappa = 2;

		//For output formatting
		this.KAPPA[0] = this.KAPPA[1] = this.kappa;

		return;
	}

	/************************************************
	* Function: Run
	* Input Parameter: seq1, seq2
	* Output: Main function for calculating Ka&Ks.
	* Return Value: void
	*************************************************/
	Run(seq1, seq2) {

		this.preProcess(seq1, seq2);

		this.S = this.L[2] / 3 + this.L[4];
		this.N = this.L[0] + 2 * this.L[2] / 3;

		this.Sd = this.L[2] * this.A[2] + this.L[4] * this.K[4];
		this.Ks = this.Sd / this.S;

		this.Nd = this.L[0] * this.K[0] + this.L[2] * this.B[2];
		this.Ka = this.Nd / this.N;

		this.t = (this.S * this.Ks + this.N * this.Ka) / (this.S + this.N);

		return this.parseOutput();
	}
}


class MLWL85 extends LWL85 {
	constructor() {
		super();
		this.name = "MLWL";
	}


	/*For more detail see reference: Tzeng Y-H, Pan R, Li W-H  (2004)  Mol. Biol. Evol.*/
	TransitionTransversion(codon1, codon2, pos) {

		//1:synonymous, 0:nonsynonymous, -1:uncalculate
		let isSyn = -1;

		//Ile: ATT, ATC, ATA; Met: ATA
		if ((codon1 == "ATA" && codon2 == "ATG" && pos == 2) || (codon1 == "ATG" && codon2 == "ATA" && pos == 2)) {
			isSyn = 0;
		}
		if ((codon1 == "ATA" && (codon2 == "ATC" || codon2 == "ATT") && pos == 2) || ((codon1 == "ATC" || codon1 == "ATT") && codon2 == "ATA" && pos == 2)) {
			isSyn = 1;
		}

		//Arg: CGA, CGG, AGA, AGG
		if ((codon1 == "CGA" && codon2 == "AGA" && pos == 0) || (codon1 == "AGA" && codon2 == "CGA" && pos == 0)) {
			isSyn = 1;
		}
		if ((codon1 == "CGG" && codon2 == "AGG" && pos == 0) || (codon1 == "AGG" && codon2 == "CGG" && pos == 0)) {
			isSyn = 1;
		}

		//Synonymous: A<->G, C<->T
		//Normal situation	
		if (isSyn == -1) {
			let sum = this.convertChar(codon1[pos]) + this.convertChar(codon2[pos]);
			if (sum == 5 || sum == 1)
				isSyn = 1;
			else
				isSyn = 0;
		}

		let class1 = this.getCodonClass(codon1, pos);
		let class2 = this.getCodonClass(codon2, pos);
		if (isSyn == 1) {
			this.Si_temp[class1] += 0.5;
			this.Si_temp[class2] += 0.5;
		}
		if (isSyn == 0) {
			this.Vi_temp[class1] += 0.5;
			this.Vi_temp[class2] += 0.5;
		}

		return 0;
	}

	/* One of differences between MLWL85 and LWL85 is allowing for kappa in S and N */
	Run(stra, strb) {

		let i;
		let ts = 0.0, tv = 0.0;	//Transition, Transversion
		let ai = new Array(5).fill(0), bi = new Array(5).fill(0);

		for (i = 0; i < stra.length; i += 3) {
			this.CountSiteAndDiff(stra.substr(i, 3), strb.substr(i, 3));
		}

		for (i = 0; i < 5; i += 2) {

			ai[i] = bi[i] = 0;

			ts += this.Si[i];
			tv += this.Vi[i];

			this.L[i] = this.L[i] / 2.0;
			this.Pi[i] = this.Si[i] / this.L[i];
			this.Qi[i] = this.Vi[i] / this.L[i];


			if ((1 - 2 * this.Pi[i] - this.Qi[i]) > 0 && (1 - 2 * this.Qi[i]) > 0) {
				ai[i] = 1 / (1 - 2 * this.Pi[i] - this.Qi[i]);
				if (this.GAMMA == 4 || (Math.abs(this.GAMMA - 0.6) < SMALLVALUE) | (this.GAMMA == -1)) {
					this.name = "GMLWL";
				}
				if (this.GAMMA == 4 || (Math.abs(this.GAMMA - 0.6) < SMALLVALUE)) {
					ai[i] = 1 - 2 * this.Pi[i] - this.Qi[i];
					ai[i] = Math.pow(ai[i], -1.0 / this.GAMMA) - 1;
				}
				bi[i] = 1 / (1 - 2 * this.Qi[i]);
				if (this.GAMMA == 4 || (Math.abs(this.GAMMA - 0.6) < SMALLVALUE)) {
					bi[i] = 1 - 2 * this.Qi[i];
					bi[i] = Math.pow(bi[i], -1.0 / this.GAMMA) - 1;

				}

				//zhangyubin add 
			}


			if (this.GAMMA == 4 || (Math.abs(this.GAMMA - 0.6) < SMALLVALUE)) {
				this.B[i] = this.GAMMA * 0.5 * bi[i];
				this.A[i] = this.GAMMA * 0.5 * ai[i] - 0.25 * this.GAMMA * bi[i];

				this.K[i] = this.A[i] + this.B[i];


			}
			else {
				if (ai[i] > 0 && bi[i] > 0) {

					if (Math.log(bi[i]) >= 0) {
						this.B[i] = 0.5 * Math.log(bi[i]);
					}

					if ((0.5 * Math.log(ai[i]) - 0.25 * Math.log(bi[i])) >= 0) {
						this.A[i] = 0.5 * Math.log(ai[i]) - 0.25 * Math.log(bi[i]);
					}

					this.K[i] = this.A[i] + this.B[i];
				}

			}
		}

		this.kappa = 2.0 * ts / tv;
		if (ts < SMALLVALUE || tv < SMALLVALUE) this.kappa = 1;

		this.KAPPA[0] = this.KAPPA[1] = this.kappa;

		if (this.kappa > 2.0) {
			this.S = (this.kappa - 1) * this.L[2] / (this.kappa + 1) + this.L[4];
			this.N = this.L[0] + 2 * this.L[2] / (this.kappa + 1);
		}
		else {
			if (this.kappa > 0.5) {
				this.S = (this.kappa - 0.5) * this.L[2] / (this.kappa + 1.5) + this.L[4];
				this.N = this.L[0] + 2 * this.L[2] / (this.kappa + 1.5);
			}
			else {
				this.S = this.L[2] / 3 + this.L[4];
				this.N = 2 * this.L[2] / 3 + this.L[0];
			}
		}

		this.Sd = this.L[2] * this.A[2] + this.L[4] * this.K[4];
		this.Ks = this.Sd / this.S;

		this.Nd = this.L[0] * this.K[0] + this.L[2] * this.B[2];
		this.Ka = this.Nd / this.N;

		this.t = (this.S * this.Ks + this.N * this.Ka) / (this.S + this.N);

		return this.parseOutput();
	}
}
