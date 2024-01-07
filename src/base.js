/************************************************************
* Filename: base.js
* Abstract: Definition of base class for KaKs methods.

* Version: js-3.1
* Author: Fanix
* Date: December.26, 2023

*************************************************************/

function gammap(x, alpha) { return alpha * (1 - x ** (-1 / alpha)) }
/******** Global variables ********/


/*						The Genetic Codes 
http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c
	Last update of the Genetic Codes: October 06, 2022 */
/* Genetic standard codon table, ! =stop codon */
let translation = {
	1: "FFLLSSSSYY!!CC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	2: "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS!!VVVVAAAADDEEGGGG",
	3: "FFLLSSSSYY!!CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	4: "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	5: "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
	6: "FFLLSSSSYYQQCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	9: "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
	10: "FFLLSSSSYY!!CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	11: "FFLLSSSSYY!!CC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	12: "FFLLSSSSYY!!CC!WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	13: "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
	14: "FFLLSSSSYYY!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
	15: "FFLLSSSSYY!QCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	16: "FFLLSSSSYY!LCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	21: "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
	22: "FFLLSS!SYY!LCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	23: "FF!LSSSSYY!!CC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	24: "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
	25: "FFLLSSSSYY!!CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	26: "FFLLSSSSYY!!CC!WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	27: "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	28: "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	29: "FFLLSSSSYYYYCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	30: "FFLLSSSSYYEECC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	31: "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	33: "FFLLSSSSYYY!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG"
}
//********End of Global variables**********

//Constructor function
class Base {
	constructor(genetic_code = 1) {
		this.translation = translation[genetic_code];
		this.DNA = "TCAG";
		this.encodedDNA = [0, 1, 2, 3];
		this.DNAencode = { 'T': 0, 'U': 0, 'C': 1, 'A': 2, 'G': 3 };
	}

	/**
	* encode DNA strings (alphabet to number)
	* @param {string} seq DNA sequence to be encoded
	* @returns {string} encoded sequence
	*/
	encodeDNA(seq) {
		return seq.split("").map(x => this.DNAencode[x] ?? "N").join("");
	}

	/**
	* decode DNA strings (number to alphabet)
	* @param {string} seq encoded DNA sequence
	* @returns {string} deencod sequence
	*/
	decodeDNA(seq) {
		return seq.split("").map(x => this.DNA[x]).join("");
	}

	/**
	* encode codon (triplet to number)
	* @param {string} codon encoded DNA triplet
	* @returns {number} codon ID
	*/
	encodeCodon(codon) {
		// if(typeof codon==="number") return codon;
		return parseInt(codon, 4);
	}

	/**
	* decode codon (number to triplet)
	* @param {number} codonID codon ID
	* @returns {string} codon (encoded DNA triplet)
	*/
	decodeCodon(codonID) {
		return (+codonID).toString(4).padStart(3, 0);
	}

	/**
	* convert (encoded) codon to amino acid
	* @param {string|number} codon encoded DNA triplet or codon ID
	* @returns {string} codon ID
	*/
	getAminoAcid(codon) {
		if (typeof codon === "number") return this.translation[codon];
		if (typeof codon === "string") return this.translation[this.encodeCodon(codon)];
		return
	}

	/**
	* test if two condon is synonyms
	* @param {string|number} codon1
	* @param {string|number} codon2
	* @returns {boolean}
	*/
	isSynonyms(codon1, codon2) {
		return this.getAminoAcid(codon1) == this.getAminoAcid(codon2);
	}

	/**
	* apply a mutation on a sequence at a specific position
	* @param {string} seq sequence to be mutated
	* @param {number} pos position to be mutated
	* @param {char} replacement mutate this position to
	* @returns {string} mutated string
	*/
	mutate(seq, pos, replacement) {
		return seq.slice(0, pos) + replacement + seq.slice(pos + 1);
	}

	/**
	* enumerate all possible mutants of a (encoded) sequence at a specific position
	* @param {string} seq sequence to be mutated
	* @param {number} pos position to be mutated
	* @returns {string[]} mutated string
	*/
	getMutatedCodonList(seq, pos) {
		return this.encodedDNA
			.filter(x => x != seq[pos])
			.map(x => this.mutate(seq, pos, x));
	}

	/**
	* get transition base
	* @param {string} nucl encoded base
	* @returns {string} encoded transition base
	*/
	transite(nucl) {
		return {
			"2": "3", "3": "2", "1": "0", "0": "1",
		}[nucl];
	}

	/**
	* get complementary base
	* @param {string} nucl encoded base
	* @returns {string} encoded complementary base
	*/
	complement(nucl) {
		return {
			"2": "0", "0": "2", "1": "3", "3": "1",
		}[nucl];
	}
}

MatrixUtils = {
	/**
	* create identity matrix
	* @param {int} n dimension of them matrix
	* @returns {int[][]} a n*n identity matrix of
	*/
	diag: function (n) {
		let A = new Array(n).fill(0).map(_ => new Array(n).fill(0));
		for (let i = 0; i < n; i++) A[i][i] = 1;
		return A;
	},
}

VectorUtils = {
	/**
	* Sum of vector elements
	* @param {number[]} v vector (numeric array)
	* @returns {number} sum of vector elements
	*/
	sum: function (v) {
		return v.reduce((a, x) => a + x);
	},

	/**
	* Norm of a vector
	* @param {number[]} v vector (numeric array)
	* @returns {number} the 2-norm of the vector
	*/
	norm: function (v) {
		return Math.sqrt(v.reduce((a, x) => a + x ** 2, 0));
	},

	/**
	* Inner product of two vectors
	* @param {number[]} v vector (numeric array)
	* @param {number[]} u vector (numeric array)
	* @returns {number} the inner product of `v` and `u`
	*/
	inner: function (v, u) {
		let s = 0;
		for (let i in v) s += v[i] * u[i];
		return s;
	},

	/**
	* Vector sum
	* @param {number[]} v vector (numeric array)
	* @param {number[]} u vector (numeric array)
	* @returns {number[]} the vectorial sum of `v` and `u`
	*/
	add: function (v, u) {
		return v.map((x, i) => x + u[i]);
	},

	/**
	* Vector difference
	* @param {number[]} v vector (numeric array)
	* @param {number[]} u vector (numeric array)
	* @returns {number[]} the vectorial difference of `v` and `u`
	*/
	sub: function (v, u) {
		return v.map((x, i) => x - u[i]);
	},

	/**
	* Vector multiplication
	* @param {number[]} v vector (numeric array)
	* @param {number} a multiplier
	* @returns {number[]} `v` multiplied by `a`
	*/
	mult: function (v, a) {
		return v.map(x => a * x);
	},

	/**
	* Integer sequence generator
	* @param {int} n sequence length
	* @param {int} init initial value
	* @returns {int[]}
	*/
	seq: function (n, init = 0) {
		return new Array(n).fill(0).map((_, i) => i + init);
	}
}
