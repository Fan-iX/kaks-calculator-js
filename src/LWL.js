/*******************************************************************
* Filename: LWL.js
* Abstract: Definition of LWL and LPB class.

* Version: js-3.1
* Author: Fanix
* Date: December.22, 2023

* Re-implementation KaKs_Calculator/LWL85 and KaKs_Calculator/LPB93

  References: 
  Li WH, Wu CI, Luo CC  (1985)  A new method for
  estimating synonymous and nonsynonymous rates of nucleotide 
  substitution considering the relative likelihood of nucleotide
  and codon changes. Mol. Biol. Evol. 2:150-174.

  Li WH  (1993)  Unbiased estimation of the Rates of synonymous
  and nonsynonymous substitution. J. Mol. Evol. 36:96-99.

  Pamilo P, Bianchi NO  (1993)  Evolution of the Zfx and Zfy 
  genes: rates and interdependence between the genes. Mol. Biol.
  Evol. 10:271-281.

  Tzeng Y-H, Pan R, Li W-H  (2004)  Comparison of Three Methods
  for Estimating Rates of Synonymous and Nonsynonymous Nucleotide
  Substitutions. Mol. Biol. Evol. 21:2290-2298.
******************************************************************/

class LWL extends Base {
	constructor(genetic_code = 1) {
		super(genetic_code)
		this.method = "LWL";
	}

	/**
	* test if two nucleotide is transition
	* @param {string} nucl1 encoded nucleotide
	* @param {string} nucl2 encoded nucleotide
	* @returns {boolean}
	*/
	isTransition(nucl1, nucl2) {
		return nucl1 == this.transite(nucl2)
	}

	/**
	* Count 0,2,4-fold sites
	* @param {string} seq encoded DNA sequence
	* @returns {{0: number, 2: number, 4: number}} number of 0,2,4-fold sites
	*/
	CountSites(seq) {
		let L = { 0: 0, 2: 0, 4: 0 }
		for (let codon of seq.match(/.{3}/g))
			for (let i of [0, 1, 2])
				L[this.getCodonClass(codon, i)]++
		return L
	}

	/**
	* Count synonymous(Si) and nonsynonymous(Vi) substitutions
	* @param {string} seq1 encoded DNA sequence
	* @param {string} seq2 encoded DNA sequence
	* @returns {[{0: number, 2: number, 4: number}, {0: number, 2: number, 4: number}]} `[Si, Vi]` number of 0,2,4-fold synonymous and nonsynonymous substitutions
	*/
	CountDiff(seq1, seq2) {
		let arr1 = seq1.match(/.{3}/g),
			arr2 = seq2.match(/.{3}/g),
			Si = { 0: 0, 2: 0, 4: 0 },
			Vi = { 0: 0, 2: 0, 4: 0 }

		for (let h in arr1) {
			let nstop = 0;
			let codon1 = arr1[h], codon2 = arr2[h];

			//Count differences
			let diff = [0, 1, 2].filter(i => codon1[i] != codon2[i]),
				ndiff = diff.length,
				//the number of evolution pathway.
				npath = { 1: 1, 2: 2, 3: 6 }[ndiff];

			//Two codons are identical
			if (ndiff == 0) continue;

			this.Si_temp = { 0: 0, 2: 0, 4: 0 };
			this.Vi_temp = { 0: 0, 2: 0, 4: 0 };

			//One difference
			if (ndiff == 1) {
				this.TransitionTransversion(codon1, codon2, diff[0]);
			}

			//Two differences: 2 pathways for evolution (i,j)
			if (ndiff == 2) {
				for (let i of diff) {
					for (let j of diff) {
						if (i != j) {
							let c = this.mutate(codon1, i, codon2[i]);
							if (this.getAminoAcid(c) != '!') {
								//pathway: codon1 <-> temp1 <-> codon2
								this.TransitionTransversion(codon1, c, i);
								this.TransitionTransversion(c, codon2, j);
							}
							else nstop++;
						}
					}
				}
			}

			//Three differences: 6 pathways for evolution (i,j,k)
			if (ndiff == 3) {
				for (let i of diff) {
					for (let j of diff) {
						for (let k of diff) {
							if ((i != j) && (i != k) && (j != k)) {
								let c1 = this.mutate(codon1, i, codon2[i]),
									c2 = this.mutate(c1, j, codon2[j]);
								if (this.getAminoAcid(c1) != '!' && this.getAminoAcid(c2) != '!') {
									//pathway: codon1 <-> temp1 <-> temp2 <-> codon2
									this.TransitionTransversion(codon1, c1, i);
									this.TransitionTransversion(c1, c2, j);
									this.TransitionTransversion(c2, codon2, k);
								}
								else nstop++;
							}
						}
					}
				}
			}

			//Add pair-codon's differences to Si and Vi
			if (npath - nstop > 0)
				for (let i of [0, 2, 4]) {
					Si[i] += (this.Si_temp[i] / (npath - nstop));
					Vi[i] += (this.Vi_temp[i] / (npath - nstop));
				}
		}
		return [Si, Vi]
	}

	/**
	* Fold class of codon substitution
	* @param {string} codon1 encoded codon
	* @param {string} codon2 encoded codon
	* @returns {0|2|4} class of codon substitution
	*/
	getCodonClass(codon, pos) {
		let codonClass = this.getMutatedCodonList(codon, pos).
			filter(c => this.getAminoAcid(c) != '!' && this.isSynonyms(c, codon)).length
		if (codonClass == 1) codonClass = 2;
		else if (codonClass == 3) codonClass = 4;

		return codonClass;
	}

	/**
	* Sum up synonymous and nonsynonymous substitutions of two codons at a given position.
	* @param {string} codon1 encoded codon
	* @param {string} codon2 encoded codon
	*/
	TransitionTransversion(codon1, codon2, pos) {
		//1:synonymous, 0:nonsynonymous
		let class1 = this.getCodonClass(codon1, pos);
		let class2 = this.getCodonClass(codon2, pos);

		//Arg: CGA(132), CGG(133), AGA(232), AGG (233)
		if (
			(codon1 == "132" && codon2 == "232" && pos == 0) ||
			(codon1 == "232" && codon2 == "132" && pos == 0) ||
			(codon1 == "133" && codon2 == "233" && pos == 0) ||
			(codon1 == "233" && codon2 == "133" && pos == 0)
		) {
			this.Si_temp[class1] += 0.5;	//different from the following
			this.Vi_temp[class2] += 0.5;
			return
		}

		/* Normal situation */
		//Synonymous: T(0)<->C(1), A(2)<->G(3) 	
		if (this.isTransition(codon1[pos], codon2[pos])) {
			this.Si_temp[class1] += 0.5;
			this.Si_temp[class2] += 0.5;
		} else {
			this.Vi_temp[class1] += 0.5;
			this.Vi_temp[class2] += 0.5;
		}
	}

	/**
	* Correct multiple substitutions
	* @param {number} P synonymous substitutions rate
	* @param {number} Q nonsynonymous substitutions rate
	* @returns {[number, number]} corrected substitutions rate
	*/
	CorrectFreq(P, Q) {
		let ai = 1 - 2 * P - Q, bi = 1 - 2 * Q;
		ai = -Math.log(ai); bi = -Math.log(bi);
		return [0.5 * ai - 0.25 * bi, 0.5 * bi]
	}

	/**
	* calculate Ka/Ks
	* @param {string} seq1 encoded DNA sequence1
	* @param {string} seq2 encoded DNA sequence2
	* @returns {{Ka: number, Ks: number, S: number, N: number, Sd: number, Nd: number, t: number, method: string, KAPPA: number[], Si:{0, 2, 4}, Vi:{0, 2, 4}, L:{0, 2, 4}}}
	*/
	calc(seq1, seq2) {
		let A = {}, B = {}
		let L = this.CountSites(seq1 + seq2);
		for (let i of [0, 2, 4]) L[i] /= 2;
		let [Si, Vi] = this.CountDiff(seq1, seq2);
		let ts = VectorUtils.sum(Object.values(Si)),
			tv = VectorUtils.sum(Object.values(Vi));

		for (let i of [0, 2, 4])
			[A[i], B[i]] = this.CorrectFreq(Si[i] / L[i], Vi[i] / L[i]);
		let kappa = tv > 0 ? 2 * ts / tv : 2

		let S = L[2] / 3 + L[4],
			N = L[0] + 2 * L[2] / 3,
			Sd = L[2] * A[2] + L[4] * (A[4] + B[4]),
			Nd = L[0] * (A[0] + B[0]) + L[2] * B[2],
			Ks = Sd / S,
			Ka = Nd / N
		return {
			Ka, Ks, S, N, Sd, Nd,
			t: (S * Ks + N * Ka) / (S + N),
			KAPPA: [kappa, kappa, 1, 1, 1, 1],
			Si, Vi, L,
			method: this.method
		};
	}
}

class MLWL extends LWL {
	constructor(genetic_code = 1) {
		super(genetic_code);
		this.method = "MLWL";
	}

	/**
	* Sum up synonymous and nonsynonymous substitutions of two codons at a given position, see *Tzeng Y-H, et al, 2004*
	* @param {string} codon1 encoded codon
	* @param {string} codon2 encoded codon
	*/
	TransitionTransversion(codon1, codon2, pos) {
		//1:synonymous, 0:nonsynonymous, -1:uncalculate
		let isSyn = -1;

		//Ile: ATT(200), ATC(201), ATA(202); Met: ATA(203)
		if (
			(codon1 == "202" && codon2 == "203" && pos == 2) ||
			(codon1 == "203" && codon2 == "202" && pos == 2)
		) {
			isSyn = 0;
		}
		if ((codon1 == "202" && (codon2 == "201" || codon2 == "200") && pos == 2) ||
			((codon1 == "201" || codon1 == "200") && codon2 == "202" && pos == 2)) {
			isSyn = 1;
		}

		//Arg: CGA(132), CGG(133), AGA(232), AGG (233)
		if (
			(codon1 == "132" && codon2 == "232" && pos == 0) ||
			(codon1 == "232" && codon2 == "132" && pos == 0) ||
			(codon1 == "133" && codon2 == "233" && pos == 0) ||
			(codon1 == "233" && codon2 == "133" && pos == 0)
		) {
			isSyn = 1;
		}

		//Synonymous: A<->G, C<->T
		//Normal situation	
		if (isSyn == -1) {
			if (this.isTransition(codon1[pos], codon2[pos]))
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
	}

	/**
	* calculate Ka/Ks
	* @param {string} seq1 encoded DNA sequence1
	* @param {string} seq2 encoded DNA sequence2
	* @returns {{Ka: number, Ks: number, S: number, N: number, Sd: number, Nd: number, t: number, method: string, KAPPA: number[], Si:{0, 2, 4}, Vi:{0, 2, 4}, L:{0, 2, 4}}}
	*/
	calc(seq1, seq2) {
		let A = {}, B = {}

		let L = this.CountSites(seq1 + seq2);
		for (let i of [0, 2, 4]) L[i] /= 2;
		let [Si, Vi] = this.CountDiff(seq1, seq2);
		let ts = VectorUtils.sum(Object.values(Si)),
			tv = VectorUtils.sum(Object.values(Vi));

		for (let i of [0, 2, 4])
			[A[i], B[i]] = this.CorrectFreq(Si[i] / L[i], Vi[i] / L[i]);

		let kappa = 2.0 * ts / tv;
		if (ts < 1e-6 || tv < 1e-6) kappa = 1;

		let S, N
		if (kappa > 2.0) {
			S = (kappa - 1) * L[2] / (kappa + 1) + L[4];
			N = L[0] + 2 * L[2] / (kappa + 1);
		}
		else {
			if (kappa > 0.5) {
				S = (kappa - 0.5) * L[2] / (kappa + 1.5) + L[4];
				N = L[0] + 2 * L[2] / (kappa + 1.5);
			}
			else {
				S = L[2] / 3 + L[4];
				N = 2 * L[2] / 3 + L[0];
			}
		}

		let Sd = L[2] * A[2] + L[4] * (A[4] + B[4]),
			Ks = Sd / S,
			Nd = L[0] * (A[0] + B[0]) + L[2] * B[2],
			Ka = Nd / N,
			t = (Sd + Nd) / (S + N)

		return {
			Ka, Ks, S, N, Sd, Nd, t,
			KAPPA: [kappa, kappa, 1, 1, 1, 1],
			Si, Vi, L,
			method: this.method
		};
	}
}

class LPB extends LWL {
	constructor(genetic_code = 1) {
		super(genetic_code);
		this.method = "LPB";
		this.Si_temp = { 0: 0, 2: 0, 4: 0 };
		this.Vi_temp = { 0: 0, 2: 0, 4: 0 };
	}

	/**
	* calculate Ka/Ks
	* @param {string} seq1 encoded DNA sequence1
	* @param {string} seq2 encoded DNA sequence2
	* @returns {{Ka: number, Ks: number, S: number, N: number, Sd: number, Nd: number, t: number, method: string, KAPPA: number[], Si:{0, 2, 4}, Vi:{0, 2, 4}, L:{0, 2, 4}}}
	*/
	calc(seq1, seq2) {
		let A = {}, B = {}, K = {}
		let L = this.CountSites(seq1 + seq2);
		for (let i of [0, 2, 4]) L[i] /= 2;
		let [Si, Vi] = this.CountDiff(seq1, seq2);
		let ts = VectorUtils.sum(Object.values(Si)),
			tv = VectorUtils.sum(Object.values(Vi));

		for (let i of [0, 2, 4]) {
			[A[i], B[i]] = this.CorrectFreq(Si[i] / L[i], Vi[i] / L[i]);
			K[i] = A[i] + B[i]
		}

		//For output formatting
		let kappa = tv > 0 ? 2 * ts / tv : 2;

		let Ks = B[4] + (L[2] * A[2] + L[4] * A[4]) / (L[2] + L[4]),
			Ka = A[0] + (L[0] * B[0] + L[2] * B[2]) / (L[0] + L[2]),
			Sd = L[2] * A[2] + L[4] * K[4],
			Nd = L[2] * B[2] + L[0] * K[0],
			S = Sd / Ks,
			N = Nd / Ka,

			t = (L[0] * K[0] + L[2] * K[2] + L[4] * K[4]) / (L[0] + L[2] + L[4]);

		return {
			Ka, Ks, S, N, Sd, Nd, t,
			KAPPA: [kappa, kappa, 1, 1, 1, 1],
			Si, Vi, L,
			method: this.method
		};
	}
}

class MLPB extends LPB {
	constructor(genetic_code = 1) {
		super(genetic_code);
		this.method = "MLPB";
	}

	/**
	* Sum up synonymous and nonsynonymous substitutions of two codons at a given position, see *Tzeng Y-H, et al, 2004*
	* @param {string} codon1 encoded codon
	* @param {string} codon2 encoded codon
	*/
	TransitionTransversion(codon1, codon2, pos) {
		//1:synonymous, 0:nonsynonymous, -1:uncalculate
		let isSyn = -1;

		//Ile: ATT(200), ATC(201), ATA(202); Met: ATA(203)
		if (
			(codon1 == "202" && codon2 == "203" && pos == 2) ||
			(codon1 == "203" && codon2 == "202" && pos == 2)
		) {
			isSyn = 0;
		}
		if ((codon1 == "202" && (codon2 == "201" || codon2 == "200") && pos == 2) ||
			((codon1 == "201" || codon1 == "200") && codon2 == "202" && pos == 2)) {
			isSyn = 1;
		}

		//Arg: CGA(132), CGG(133), AGA(232), AGG (233)
		if (
			(codon1 == "132" && codon2 == "232" && pos == 0) ||
			(codon1 == "232" && codon2 == "132" && pos == 0) ||
			(codon1 == "133" && codon2 == "233" && pos == 0) ||
			(codon1 == "233" && codon2 == "133" && pos == 0)
		) {
			isSyn = 1;
		}

		//Synonymous: A<->G, C<->T
		//Normal situation	
		if (isSyn == -1) {
			if (this.isTransition(codon1[pos], codon2[pos]))
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
	}
}
