/************************************************************
* Filename: MA.js
* Abstract: Definition of model-averged methods' classes.

* Version: js-3.1
* Author: Fanix
* Date: January.6, 2024

* Re-implementation of KaKs_Calculator/MSMA

  References:
  Goldman N, Yang Z  (1994)  A codon-based model of nucleotide
  substitution for protein-coding DNA sequences. Mol. Biol. 
  Evol. 11:725-736.

  Posada, D. and Buckley, T.R. (2004) Model Selection and Model Averaging
  in Phylogenetics: Advantages of Akaike Information Criterion and Bayesian
  Approaches over Likelihood Ratio Tests, Syst. Biol., 53, 793-808.

  Sullivan, J. and Joyce, P. (2005) Model Selection in Phylogenetics, 
  Annual Review of Ecology, Evolution, and Systematics, 36, 445-466.  
*************************************************************/

class MA extends GY {
	constructor(genetic_code = 1) {
		super(genetic_code = 1);
		this.genetic_code = genetic_code
		this.method = "MA";
	}

	calc(seq1, seq2) {
		let seq_length = seq1.length;
		let snp = VectorUtils.seq(seq_length).filter(i => seq1[i] != seq2[i]).length;
		let fp = {}
		let arr1 = seq1.match(/.{3}/g),
			arr2 = seq2.match(/.{3}/g)

		for (let h in arr1) {
			let codon1 = arr1[h], codon2 = arr2[h]
			if (codon1 < codon2) [codon1, codon2] = [codon2, codon1]
			if (!fp[codon1]) fp[codon1] = {}
			if (!fp[codon1][codon2]) fp[codon1][codon2] = 0
			fp[codon1][codon2]++
		}

		let result = ["JC", "F81", "K2P", "HKY", "TNEF", "TN", "K3P", "K3PUF", "TIMEF", "TIM", "TVMEF", "TVM", "SYM", "GTR"].map(model => new GY(this.genetic_code, model).calc(seq1, seq2))
		let minAICc = Math.min(...result.map(x => x.AICc)),
			diff = result.map(x => x.AICc - minAICc),
			w = diff.map(x => 1 / VectorUtils.sum(diff.map(y => Math.exp(0.5 * x - 0.5 * y))))
		let I = [
			[0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0],
			[1, 1, 0, 0, 0, 0], [1, 1, 0, 0, 0, 0],
			[1, 1, 0, 0, 0, 0], [1, 1, 0, 0, 0, 0],
			[1, 1, 1, 1, 0, 0], [1, 1, 1, 1, 0, 0],
			[1, 1, 1, 1, 0, 0], [1, 1, 1, 1, 0, 0],
			[1, 1, 1, 1, 1, 0], [1, 1, 1, 1, 1, 0],
			[1, 1, 1, 1, 1, 0], [1, 1, 1, 1, 1, 0],
		]
		let t = VectorUtils.sum(w.map((_, i) => w[i] * result[i].t * 3))
		let omega = VectorUtils.sum(w.map((_, i) => w[i] * result[i].Ka / result[i].Ks))
		let KAPPA = new Array(6).fill(0),
			sum = new Array(6).fill(0);
		result.forEach((res, i) => res.KAPPA.forEach((kap, j) => {
			KAPPA[j] += w[i] * I[i][j] * kap;
			sum[j] += w[i] * I[i][j];
		}))
		KAPPA = KAPPA.map((_, i) => sum[i] == 0 ? 1 : KAPPA[i] / sum[i])

		let efpi = new Array(64).fill(0);
		this.FROM61.forEach(i => efpi[i] = 1 / this.FROM61.length);
		let ufpi = this.getFreqency(seq1 + seq2)
		let pi = new Array(64).fill(0);
		w.forEach((_, i) => (i % 2 == 0 ? efpi : ufpi).forEach((p, j) => pi[j] += w[i] * p))


		let [rs, ra, rs0, ra0] = this.EigenQc(omega, KAPPA, pi);

		let S = rs0 / (rs0 + ra0) * seq_length, N = seq_length - S,
			Sd = snp * rs / (rs + ra), Nd = snp - Sd,
			Ks = t / 3 * rs / (rs + ra) * (rs0 + ra0) / rs0,
			Ka = t / 3 * ra / (rs + ra) * (rs0 + ra0) / ra0

		let n = this.FROM61.length;
		let [Root, U, V] = this.GetEigenQ(omega, KAPPA, pi);
		let expt = Root.map(r => Math.exp(t * r));
		let lnL = 0
		for (let c1 in fp) {
			let z0 = this.FROM64[this.encodeCodon(c1)]
			for (let c2 in fp[c1]) {
				let z1 = this.FROM64[this.encodeCodon(c2)]
				let fh = 0;
				for (let k = 0; k < n; k++) {
					fh += U[z0 * n + k] * expt[k] * V[k * n + z1];
				}
				fh *= pi[this.encodeCodon(c1)];
				lnL += Math.log(fh) * fp[c1][c2];
			}
		}

		return {
			Ka, Ks, S, N, Sd, Nd,
			t: (S * Ks + N * Ka) / (S + N),
			lnL, KAPPA,
			method: this.method
		};
	}
}
