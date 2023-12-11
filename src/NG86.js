/*********************************************************
* Filename: NG86.js
* Abstract: Definition of NG86 class.

* Version: js-2.0
* Author: Fanix
* Date: December.11, 2022

* Adapted from: KaKs_Calculator/NG86.cpp and KaKs_Calculator/NG86.h
* Copyright (C) 2009, Beijing Institute of Genomics of CAS
* All rights reserved.

* Version: 2.0
* Author: Da-Peng Wang(wangdp@big.ac.cn), Yu-Bin Zhang (ybzhang@big.ac.cn)
* Date: Jun.1, 2009

* Version: 1.0
* Author: Zhang Zhang (zhang.zhang@yale.edu)
* Date: Feb.21, 2005

  Reference: Nei M, Gojobori T  (1986)  Simple methods for
  estimating the numbers of synonymous and nonsynonymous 
  nucleotide substitutions. Mol Biol Evol 3:418-426.
**********************************************************/

class NG86 extends Base {
	constructor() {
		super();
		this.name = "NG";
		this.Ps = 0;
		this.Pn = 0;
	}

	getCondonSite(codon) {
		let i, j, stop;
		let temp = "";
		let syn = 0.0;

		if (this.getAminoAcid(codon) == '!') return;

		/* Synonymous sites only occur at first and third position in a codon */
		for (i = 0, stop = 0; i < 3; i += 2) {
			for (j = 0; j < 4; j++) {
				temp = codon.split("");
				if (j != this.convertChar(temp[i])) {
					temp[i] = this.convertInt(j);
					if (this.getAminoAcid(temp.join("")) == '!') {
						stop++;
					}
					else {
						if (this.getAminoAcid(temp.join("")) == this.getAminoAcid(codon))
							syn++;
					}
				}
			}
		}
		this.S += (syn / 3.0);
		this.N += (3 - stop / 3.0 - syn / 3.0);
	}

	getCondonDifference(codon1, codon2) {

		let i, j, k, diff = new Array(CODONLENGTH);
		let num = 0;
		let stop = 0;
		let path = 1;
		let sd_temp = 0;
		let nd_temp = 0;
		let temp1, temp2;


		if (this.getAminoAcid(codon1) == '!' || this.getAminoAcid(codon2) == '!')
			return;

		for (i = 0; i < CODONLENGTH; i++) {
			diff[i] = -1;
			if (codon1[i] != codon2[i])
				diff[num++] = i;
		}

		//two codons are same
		if (num == 0) return;

		this.snp += num;

		//Pathway of evolution from the differences
		for (i = 1; i <= num; i++)
			path *= i;

		//Only one difference between two codons
		if (num == 1) {
			if (this.getAminoAcid(codon1) == this.getAminoAcid(codon2))
				sd_temp++;
			else
				nd_temp++;
		}

		//Two differences between two codons
		if (num == 2) {
			for (i = 0; i < num; i++)
				for (j = 0; j < num; j++)
					if (i != j) {
						temp1 = codon1.split("");
						temp1[diff[i]] = codon2[diff[i]];
						temp1 = temp1.join("")
						if (this.getAminoAcid(temp1) != '!') {

							//codon1<->temp1
							if (this.getAminoAcid(temp1) == this.getAminoAcid(codon1)) sd_temp++;
							else nd_temp++;

							//temp1<->codon2
							if (this.getAminoAcid(temp1) == this.getAminoAcid(codon2)) sd_temp++;
							else nd_temp++;
						}
						else
							stop++;
					}
		}
		//Three differences between two codons
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
								//codon1<->temp1
								if (this.getAminoAcid(temp1) == this.getAminoAcid(codon1)) sd_temp++;
								else nd_temp++;

								//temp1<->temp2
								if (this.getAminoAcid(temp2) == this.getAminoAcid(temp1)) sd_temp++;
								else nd_temp++;

								//temp2<->codon2
								if (this.getAminoAcid(codon2) == this.getAminoAcid(temp2)) sd_temp++;
								else nd_temp++;

							}
							else
								stop++;
						}
					}
				}
			}
		}
		if (path == stop) {
			//All pathways are through stop codons
			if (num == 2) {
				this.Sd += 0.5; this.Nd += 1.5;
			}
			else {
				this.Sd += 1.0; this.Nd += 2.0;
			}
		}
		else {
			this.Sd += (sd_temp / (path - stop));
			this.Nd += (nd_temp / (path - stop));
		}
	}
	PreProcess(seq1, seq2) {
		let i;

		//Count sites and differences
		for (i = 0; i < seq1.length; i = i + 3) {
			this.getCondonSite(seq1.substr(i, 3));
			this.getCondonSite(seq2.substr(i, 3));
			this.getCondonDifference(seq1.substr(i, 3), seq2.substr(i, 3));
		}

		this.S /= 2.0;
		this.N /= 2.0;

		//Scale the sum of S+N to the length of sequence.
		let y = seq1.length / (this.S + this.N);
		this.S *= y;
		this.N *= y;
	}
	kaks_formula(p) {
		//Equation (3) in the reference of NG86
		let d = 1 - (4 * p) / 3;
		if (d < 0.0) {
			d = NA;
		}
		else {
			if (this.GAMMA == 6 || this.GAMMA == -1)   //zhangyubin added
			{
				this.name = "GNG";
			}
			if (this.GAMMA == 6) {
				d = pow(d, -1.0 / 0.6) - 1;
				if (d < 0.0) //before (d>0.0)
					d = NA;
				else
					// d = (-3.0)*d/4.0;
					d = (3 * d * 0.6) / 4.0;
			}
			else {
				d = Math.log(d);
				if (d > 0.0)
					d = NA;
				else
					d = (-3.0) * d / 4.0;
			}
		}
		return d;
	}
	Run(seq1, seq2) {

		this.PreProcess(seq1, seq2);
		this.length = seq1.length
		this.Ks = this.kaks_formula(this.Sd / this.S);
		this.Ka = this.kaks_formula(this.Nd / this.N);

		this.t = (this.S * this.Ks + this.N * this.Ka) / (this.S + this.N);

		return this.parseOutput();
	}
}


/***********************************************************
  NONE: an in-house algorithm in BGI for testing Ka and Ks.
  NONE is NG86 without correction for multiple substitutions.
************************************************************/
class NONE extends NG86 {
	constructor() {
		super();
		this.name = "NONE";
	}

	Run(seq1, seq2) {
		this.PreProcess(seq1, seq2);

		this.Ks = this.Sd / this.S;
		this.Ka = this.Nd / this.N;

		this.t = (this.S * this.Ks + this.N * this.Ka) / (this.S + this.N);

		return this.parseOutput();
	}
}
