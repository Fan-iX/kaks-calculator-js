/*********************************************************
* Filename: NG.js
* Abstract: Definition of NG class.

* Version: js-3.1
* Author: Fanix
* Date: December.20, 2023

* Re-implementation of KaKs_Calculator/NG86

  Reference:
  Nei M, Gojobori T  (1986)  Simple methods for estimating
  the numbers of synonymous and nonsynonymous nucleotide
  substitutions. Mol Biol Evol 3:418-426.
**********************************************************/

class NG extends Base {
	constructor(genetic_code = 1) {
		super(genetic_code);
		this.method = "NG";
	}

	/**
	* Count synonymous(S) and nonsynonymous(N) sites
	* @param {string} seq encoded DNA sequence
	* @returns {[number, number]} `[S, N]` number of synonymous and nonsynonymous sites
	*/
	CountSites(seq) {
		let nstop = 0, nsyn = 0
		for (let codon of seq.match(/.{3}/g)) {
			for (let pos of [0, 2]) {
				for (let c of this.getMutatedCodonList(codon, pos)) {
					if (this.getAminoAcid(c) == '!') nstop++;
					else if (this.isSynonyms(codon, c)) nsyn++;
				}
			}
		}
		return [nsyn / 3, seq.length - nstop / 3 - nsyn / 3]
	}

	/**
	* Count synonymous(Sd) and nonsynonymous(Nd) substitutions
	* @param {string} seq1 encoded DNA sequence
	* @param {string} seq2 encoded DNA sequence
	* @returns {[number, number]} `[Sd, Nd]` number of synonymous and nonsynonymous substitutions
	*/
	CountDiffs(seq1, seq2) {
		let arr1 = seq1.match(/.{3}/g),
			arr2 = seq2.match(/.{3}/g);
		let Sd = 0, Nd = 0;
		for (let h in arr1) {
			let nstop = 0, sd_temp = 0, nd_temp = 0
			let codon1 = arr1[h], codon2 = arr2[h];
			if (this.getAminoAcid(codon1) == '!' || this.getAminoAcid(codon2) == '!' || codon1 == codon2) continue;
			let diff = [0, 1, 2].filter(i => codon1[i] != codon2[i]),
				ndiff = diff.length,
				npath = { 1: 1, 2: 2, 3: 6 }[diff.length];
			if (ndiff == 1) {
				if (this.isSynonyms(codon1, codon2)) sd_temp++;
				else nd_temp++;
			} else if (ndiff == 2) {
				for (let i of diff) {
					let temp1 = this.mutate(codon1, i, codon2[i])
					if (this.getAminoAcid(temp1) != '!') {
						//codon1<->temp1
						if (this.isSynonyms(temp1, codon1)) sd_temp++;
						else nd_temp++;

						//temp1<->codon2
						if (this.isSynonyms(temp1, codon2)) sd_temp++;
						else nd_temp++;
					}
					else nstop++;
				}
			} else if (ndiff == 3) {
				for (let i of diff) {
					for (let j of diff) {
						if ((i != j)) {
							let temp1 = this.mutate(codon1, i, codon2[i]);
							let temp2 = this.mutate(temp1, j, codon2[j]);
							if (this.getAminoAcid(temp1) != '!' && this.getAminoAcid(temp2) != '!') {
								//codon1<->temp1
								if (this.isSynonyms(temp1, codon1)) sd_temp++;
								else nd_temp++;

								//temp1<->temp2
								if (this.isSynonyms(temp2, temp1)) sd_temp++;
								else nd_temp++;

								//temp2<->codon2
								if (this.isSynonyms(codon2, temp2)) sd_temp++;
								else nd_temp++;

							}
							else nstop++;
						}
					}
				}
			}
			if (npath == nstop) {
				//All pathways are through stop codons
				if (ndiff == 2) {
					Sd += 0.5; Nd += 1.5;
				}
				else if (ndiff == 3) {
					Sd += 1, Nd += 2;
				}
			} else {
				Sd += sd_temp / (npath - nstop);
				Nd += nd_temp / (npath - nstop);
			}
		}
		return [Sd, Nd]
	}

	/**
	* Correct multiple substitutions, see *Nei M and Gojobori T, 1986*
	* @param {number} p substitutions rate
	* @returns {number} corrected substitutions rate
	*/
	CorrectKaksJC69(p) {
		if (p < 0) return NaN;
		return 0.75 * -Math.log(1 - p * 4 / 3);
	}

	/**
	* calculate Ka/Ks
	* @param {string} seq1 encoded DNA sequence1
	* @param {string} seq2 encoded DNA sequence2
	* @returns {{Ka: number, Ks: number, S: number, N: number, Sd: number, Nd: number, t: number, method: string}}
	*/
	calc(seq1, seq2) {
		//Count sites and differences
		let [S, N] = this.CountSites(seq1 + seq2)
		let [Sd, Nd] = this.CountDiffs(seq1, seq2)
		//Scale the sum of S+N to the length of sequence.
		let scale = seq1.length / (S + N);
		S *= scale;
		N *= scale;
		let Ks = this.CorrectKaksJC69(Sd / S),
			Ka = this.CorrectKaksJC69(Nd / N)

		return {
			Ka, Ks, S, N, Sd, Nd,
			t: (S * Ks + N * Ka) / (S + N),
			method: this.method
		};
	}
}
