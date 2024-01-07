/************************************************************
* Filename: KaKs.js
* Abstract: Declaration of KAKS class including several methods.

* Version: js-3.1
* Author: Fanix
* Date: January.6, 2024

* Adapted from: KaKs_Calculator/KaKs

*************************************************************/

class KAKS extends Base {
	constructor(genetic_code = 1) {
		super(genetic_code)
		this.genetic_code = genetic_code
		this.header = ["Sequence", "Method", "Ka", "Ks", "Ka/Ks",
			"P-Value(Fisher)", "Length", "S-Sites", "N-Sites", "Fold-Sites(0:2:4)",
			"Substitutions", "S-Substitutions", "N-Substitutions", "Fold-S-Substitutions(0:2:4)", "Fold-N-Substitutions(0:2:4)",
			"Divergence-Time", "Substitution-Rate-Ratio(rTC:rAG:rTA:rCG:rTG:rCA/rCA)", "GC(1:2:3)", "ML-Score", "AICc",
			"Akaike-Weight"];
		this.method_name = [
			"NG", "LWL", "LPB", "MLWL", "MLPB", "GY", "YN", "MYN", "GYN", "GMYN"
		];
		this.method_ref = {
			NG: "Nei, M. and Gojobori, T. (1986) Mol. Biol. Evol., 3, 418-426.",
			LWL: "Li, W.H., Wu, C.I. and Luo, C.C. (1985) Mol. Biol. Evol., 2, 150-174.",
			LPB: "Li, W.H. (1993) J. Mol. Evol., 36, 96-99.    Pamilo, P. and Bianchi, N.O. (1993) Mol. Biol. Evol., 10, 271-281.",
			MLWL: "Tzeng, Y.H., Pan, R. and Li, W.H. (2004) Mol. Biol. Evol., 21, 2290-2298.",
			MLPB: "Tzeng, Y.H., Pan, R. and Li, W.H. (2004) Mol. Biol. Evol., 21, 2290-2298.",
			GY: "Goldman, N. and Yang, Z. (1994) Mol. Biol. Evol., 11, 725-736.",
			YN: "Yang, Z. and Nielsen, R. (2000) Mol. Biol. Evol., 17, 32-43.",
			MYN: "Zhang, Z., Li, J. and Yu, J. (2006) BMC Evolutionary Biology, 6, 44.",
		};
	}

	formatResult(res) {
		let result = []
		result.push(res.seq_name, res.method)
		//Ka, Ks, Ka/Ks
		result.push(+res.Ka.toExponential(5), +res.Ks.toExponential(5), +(res.Ka / res.Ks).toExponential(5))
		//Fisher's test: p_value
		result.push(+MathUtils.fisher(res.Sd, res.Nd, res.S - res.Sd, res.N - res.Nd).toExponential(5));
		//Length of compared pairwise sequences
		result.push(res.length)
		//Synonymous(S) sites
		result.push(+res.S.toExponential(5))
		//Nonsynonymous(N) sites
		result.push(+res.N.toExponential(5))
		//L[0], L[2], L[4] only for Prof.Li's series(LWL85, LPB93...)
		result.push(res.L ? [0, 2, 4].map(x => +res.L[x].toExponential(5)).join(":") : "NA")
		//Substitutions
		result.push(res.snp)
		//Sysnonymous Substitutions(Sd)
		result.push(+res.Sd.toExponential(5))
		//Nonsysnonymous Substitutions(Nd)
		result.push(+res.Nd.toExponential(5))

		//Si, Vi for Li's series' methods(LWL85, LPB93...)
		result.push(res.Si ? [0, 2, 4].map(x => +res.Si[x].toExponential(5)).join(":") : "NA")
		result.push(res.Vi ? [0, 2, 4].map(x => +res.Vi[x].toExponential(5)).join(":") : "NA")
		//Divergence time or distance t = (S*Ks+N*Ka)/(S+N)
		result.push(+res.t.toExponential(5))
		//Substitution-Rate-Ratio(rTC:rAG:rTA:rCG:rTG:rCA/rCA)
		result.push(res.KAPPA?.map(x => +x.toExponential(5))?.join(":") || "1:1:1:1:1:1")
		//GC Content
		result.push(res.GC && +(VectorUtils.sum(res.GC) / 3).toExponential(5) + "(" + res.GC.map(x => +x.toExponential(5)).join(":") + ")")

		//Maximum Likelihood Value
		result.push(res.lnL ? +res.lnL.toExponential(5) : "NA");
		//AICc
		result.push(res.AICc ? +res.AICc.toExponential(5) : "NA");
		//Akaike weight in model selection
		result.push(res.AkaikeWeight || "NA");

		return result.join("\t") + "\n";
	}

	prepareSeq(data) {
		let [seq_name, str1, str2] = data.split('\n').map(x => x.trim());
		if (str1.length != str2.length || str1.length % 3 != 0 || str2.length % 3 != 0) {
			throw "Error. The size of two sequences in '" + seq_name + "' is not equal.";
		}
		let arr1 = this.encodeDNA(str1.toUpperCase()).match(/.{3}/g)
		let arr2 = this.encodeDNA(str2.toUpperCase()).match(/.{3}/g)
		for (let i = 0; i < arr1.length; i++) {
			if (
				/N/g.test(arr1[i]) |
				/N/g.test(arr2[i]) |
				this.getAminoAcid(arr1[i]) == "!" |
				this.getAminoAcid(arr2[i]) == "!"
			) {
				arr1[i] = arr2[i] = ""
			}
		}
		let seq1 = arr1.join(""), seq2 = arr2.join("");

		return { seq_name, seq1, seq2 }
	}

	parseAxt(axt) {
		let data = {}
		for (let seqs of axt.trim().split(/\n{2,}/)) {
			let [seq_name, seq1, seq2, remain] = seqs.split('\n').map(x => x.trim())
			if (remain) MiscUtils.error(`Error. ${seq_name} has extra lines`)
			else if (seq1.length != seq2.length) MiscUtils.error(`Error. The sequences in '${seq_name}' are not equal in length.`)
			else if (seq1.length % 3 != 0) MiscUtils.error(`Error. The sequences in '${seq_name}' are not codon-based alignment.`)
			else data[seq_name] = [seq1, seq2]
		}
		return data
	}

	Run(data, methods, args) {
		let time_start = new Date()
		console.info(`methods:${methods.join(" ")}\ngenetic_code:${this.genetic_code}`)
		let res = this.header.join("\t") + "\n";

		res += this.calculateKaKs(this.parseAxt(data), methods, args).map(x => this.formatResult(x)).join("")

		//Time used for running
		let t = (new Date() - time_start) / 1000 | 0;
		let h = t / 3600 | 0, m = (t % 3600) / 60 | 0, s = t % 60;

		//Print on display
		console.info("Mission accomplished. (Time elapsed:" + (h ? h + ":" : "") + m + ":" + s + ")");
		return res
	}

	calculateKaKs(data, methods, args) {
		let result = []
		for (let seq_name in data) {
			let [seq1, seq2] = data[seq_name].map(x => x.toUpperCase())
			let GC = this.getGCContent(seq1 + seq2)
			let arr1 = this.encodeDNA(seq1).match(/.{3}/g)
			let arr2 = this.encodeDNA(seq2).match(/.{3}/g)
			for (let i = 0; i < arr1.length; i++) {
				if (
					/N/g.test(arr1[i]) |
					/N/g.test(arr2[i]) |
					this.getAminoAcid(arr1[i]) == "!" |
					this.getAminoAcid(arr2[i]) == "!"
				) {
					arr1[i] = arr2[i] = ""
				}
			}
			seq1 = arr1.join(""); seq2 = arr2.join("");
			let length = seq1.length,
				snp = new Array(length).fill(0).filter((_, i) => seq1[i] != seq2[i]).length
			for (let method of methods) {
				let zz
				switch (method) {
					case "NG": zz = new NG(this.genetic_code); break
					case "LWL": zz = new LWL(this.genetic_code); break
					case "MLWL": zz = new MLWL(this.genetic_code); break
					case "LPB": zz = new LPB(this.genetic_code); break
					case "MLPB": zz = new MLPB(this.genetic_code); break
					case "YN": zz = new YN(this.genetic_code); break
					case "MYN": zz = new MYN(this.genetic_code); break
					case "GY":
						for (let model of args.gy_models ?? ["HKY"]) {
							zz = new GY(this.genetic_code, model);
							result.push({
								...zz.calc(seq1, seq2),
								seq_name, length, GC, snp
							})
						}
						continue
					case "MA": zz = new MA(this.genetic_code); break
					default:
						continue
				}
				result.push({
					...zz.calc(seq1, seq2),
					seq_name, length, GC, snp
				})
			}
			console.info(seq_name + " [OK]");
		}
		return result
	}

	/* Get GCC of entire sequences GC[0] and of three codon positions GC[1,2,3] */
	getGCContent(str) {
		let GCcontent = [0, 0, 0];
		str.match(/.{3}/g)
			.forEach(c => c.split("").forEach((x, i) =>
				GCcontent[i] += x == "G" || x == "C"
			))
		GCcontent = GCcontent.map(x => x / (str.length / 3))
		return GCcontent
	}
}


MathUtils = {
	/**************************************************
	* Function: fisher
	* Input Parameter: double, double, double, double
	* Output: Compute p-value by Fisher exact test
	* Return Value: double
	***************************************************/
	fisher: function (sd, nd, s, n) {
		let prob_total = 0;

		let matrix = [sd, s, nd, n],
			R = [sd + nd, s + n],
			C = [sd + s, nd + n],
			sum = sd + s + nd + n,
			numerator = 0

		//Calculate the numberator that is a constant
		numerator += this.factorial(R[0]);
		numerator += this.factorial(R[1]);
		numerator += this.factorial(C[0]);
		numerator += this.factorial(C[1]);

		//Log of Factorial of N
		let fac_sum = this.factorial(sum),
			denominator = fac_sum;
		for (let i = 0; i < 4; i++) {
			denominator += this.factorial(matrix[i]);
		}
		//Probability of current situtation
		let prob_current = Math.exp(numerator - denominator);

		//Two-tail probabilities if less than prob_current
		for (let i = 0; i <= R[0]; i++) {
			matrix[0] = i;
			matrix[1] = C[0] - i;
			matrix[2] = R[0] - i;
			matrix[3] = R[1] - C[0] + i;
			if (matrix[0] > 0 && matrix[1] > 0 && matrix[2] > 0 && matrix[3] > 0) {
				denominator = fac_sum
				for (let j = 0; j < 4; j++) {
					denominator += this.factorial(matrix[j]);
				}
				let temp = numerator - denominator;
				temp = Math.exp(numerator - denominator);
				if (temp <= prob_current) {
					prob_total += temp;
				}
			}
		}

		return prob_total;
	},
	factorial: function (n) {
		let x = 0.1659470187408462e-06 / (n + 8) +
			0.9934937113930748e-05 / (n + 7) -
			0.1385710331296526 / (n + 6) +
			12.50734324009056 / (n + 5) -
			176.6150291498386 / (n + 4) +
			771.3234287757674 / (n + 3) -
			1259.139216722289 / (n + 2) +
			676.5203681218835 / (n + 1) +
			0.9999999999995183
		return Math.log(x) - 5.58106146679532777 - (n + 1) +
			(n + 0.5) * Math.log(n + 7.5);
	}
}

MiscUtils = {
	error: function () {
		console.error(...arguments)
	}
}
