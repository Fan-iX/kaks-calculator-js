/*********************************************************
* Filename: LPB93.js
* Abstract: Definition of LPB93 and MLPB93 class.

* Version: js-2.0
* Author: Fanix
* Date: December.11, 2022

* Adapted from: KaKs_Calculator/LPB93.cpp and KaKs_Calculator/LPB93.h
* Copyright (C) 2009, Beijing Institute of Genomics of CAS
* All rights reserved.

* Version: 2.0
* Author: Da-Peng Wang(wangdp@big.ac.cn), Yu-Bin Zhang (zhangyb@big.ac.cn)
* Date: Jun.1, 2009

* Version: 1.0
* Author: Zhang Zhang (zhanghzhang@genomics.org.cn)
* Date: Jan.21, 2005

  References: 
  Li WH  (1993)  Unbiased estimation of the Rates of synonymous
  and nonsynonymous substitution. J. Mol. Evol. 36:96-99.

  Pamilo P, Bianchi NO  (1993)  Evolution of the Zfx and Zfy 
  genes: rates and interdependence between the genes. Mol. Biol.
  Evol. 10:271-281.

  Tzeng Y-H, Pan R, Li W-H  (2004)  Comparison of Three Methods
  for Estimating Rates of Synonymous and Nonsynonymous Nucleotide
  Substitutions. Mol. Biol. Evol. 21:2290-2298.
**********************************************************/


class LPB93 extends LWL85 {
	constructor() {
		super();
		this.name = "LPB";
	}

	/* Similar to LWL85 except the formulas for calculating ka and ks*/
	Run(seq1, seq2) {

		this.preProcess(seq1, seq2);

		this.Ks = this.B[4] + (this.L[2] * this.A[2] + this.L[4] * this.A[4]) / (this.L[2] + this.L[4]);
		this.Ka = this.A[0] + (this.L[0] * this.B[0] + this.L[2] * this.B[2]) / (this.L[0] + this.L[2]);

		this.Sd = this.L[2] * this.A[2] + this.L[4] * this.K[4];
		this.Nd = this.L[2] * this.B[2] + this.L[0] * this.K[0];

		this.S = this.Sd / this.Ks;
		this.N = this.Nd / this.Ka;

		this.t = (this.L[0] * this.K[0] + this.L[2] * this.K[2] + this.L[4] * this.K[4]) / (this.L[0] + this.L[2] + this.L[4]);

		return this.parseOutput();
	}
}


/*The difference between LPB93 and MLPB93 focuses on the definition of transition & transversion*/
class MLPB93 extends LPB93 {
	constructor() {
		super();
		this.name = "MLPB";
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

}
