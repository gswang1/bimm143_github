# Class 6: R Functions
Grace Wang (PID: A16968688)

- [Function basics](#function-basics)
- [Generate a DNA sequence](#generate-a-dna-sequence)
- [Generate a protein sequence](#generate-a-protein-sequence)

## Function basics

First silly function - add some numbers:

Every R function has 3 components:

- name
- input argument(s), separated by commas
- body

``` r
add <- function(x, y = 10, z = 0){
  x + y + z
}
```

Once run, I can use this function like any other function.

``` r
add(1, 100)
```

    [1] 101

``` r
add(1:4, 100)
```

    [1] 101 102 103 104

``` r
add(1)
```

    [1] 11

Functions can have “required” input arguments and “optional” input
arguments. The optional arguments are defined with `= [default]` in the
definition of the function.

``` r
add(1, 100, 10)
```

    [1] 111

## Generate a DNA sequence

> Write a function to return a DNA sequence of a user-specified length.
> Call it `generate_dna()`.

The `sample()` function can help.

``` r
bases <- c("A", "C", "G", "T")
sample(bases, size = 10, replace = T)
```

     [1] "G" "A" "A" "G" "G" "A" "A" "G" "A" "A"

Now I have a working snippet. I can put this in the body of my function.

``` r
generate_dna <- function(length = 5){
  bases <- c("A", "C", "G", "T")
  sample(bases, length, replace = T)
}

generate_dna()
```

    [1] "G" "T" "T" "C" "G"

``` r
generate_dna(100)
```

      [1] "A" "A" "G" "A" "C" "C" "T" "T" "A" "A" "T" "C" "A" "T" "G" "A" "A" "C"
     [19] "C" "G" "C" "C" "T" "A" "C" "T" "A" "G" "G" "G" "G" "T" "T" "C" "G" "T"
     [37] "A" "A" "A" "C" "A" "A" "C" "T" "G" "G" "G" "T" "T" "C" "A" "C" "C" "A"
     [55] "C" "C" "C" "T" "C" "A" "G" "T" "C" "G" "G" "A" "A" "G" "G" "T" "C" "G"
     [73] "C" "A" "C" "G" "G" "T" "C" "G" "A" "T" "T" "G" "G" "G" "G" "G" "G" "A"
     [91] "G" "C" "T" "T" "T" "C" "A" "C" "T" "G"

I want the ability to return a single sequence like “AGTACCTG” instead
of many strings in a vector.

``` r
generate_dna <- function(length = 5, together = T){
  bases <- c("A", "C", "G", "T")
  sequence <- sample(bases, length, replace = T)
  if(together){
    sequence <- paste(sequence, collapse = "")
  }
  return(sequence)
}

generate_dna()
```

    [1] "TCATC"

``` r
generate_dna(, together = F)
```

    [1] "G" "A" "A" "T" "C"

## Generate a protein sequence

> Write a function to return protein sequences of different lengths and
> test whether these sequences are unique in nature.

We can get the set of 20 canonical amino acids from the **bio3d**
package.

``` r
aa <- bio3d::aa.table$aa1[1:20]
```

> Write a protein sequence generating function that will return
> sequences of a user-specified length.

``` r
generate_protein <- function(length = 10, string = T){
  aa <- bio3d::aa.table$aa1[1:20]
  aasequence <- sample(aa, size = length, replace = T)
  if(string == T){
    aasequence <- paste(aasequence, collapse = "")
  }
  return(aasequence)
}

generate_protein(20)
```

    [1] "FRQYVLNQPGQEKDDVICSC"

> Generate random protein sequences of length 6-12 amino acids.

``` r
#generate_protein(6:12)
```

We can’t use vectors as inputs. We can fix this by editing the function
itself or by using the R **apply** family of functions.

``` r
sapply(6:12, generate_protein)
```

    [1] "HWDPMK"       "LDRMLNW"      "HPQINEHV"     "ATPLRTPKY"    "IKNAFHSAIV"  
    [6] "IWFTDTWRLHI"  "FDYNTIWSVNAV"

It would be useful to put these in FASTA format output.

``` r
seqs <- sapply(6:12, generate_protein)
seqs
```

    [1] "EDGHTN"       "EQHHYLT"      "NKRNVFDM"     "NKMPWINYT"    "GYAMSDNSLN"  
    [6] "GFCREPGKKQC"  "LMRVAYIRQLII"

``` r
ids <- paste(">ID.", 6:12, sep = "")
fasta <- paste(ids, seqs, sep = "\n")
cat(fasta, sep = "\n")
```

    >ID.6
    EDGHTN
    >ID.7
    EQHHYLT
    >ID.8
    NKRNVFDM
    >ID.9
    NKMPWINYT
    >ID.10
    GYAMSDNSLN
    >ID.11
    GFCREPGKKQC
    >ID.12
    LMRVAYIRQLII

> Determine if these sequences can be found in nature or are they
> unique? Why or why not?

This can be done using a blastp search. If both percent identity and
query cover are 100%, the sequence is not unique and is found in nature.
If the top hit has a value of \<100% for either percent identity or
query cover, then the sequence is not found in nature. For my sequences,
lengths 6 and 7 were found in nature but 8-12 were unique.

Random sequences of length \>=9 are unique and can’t be found in
databases. As such, protein sequences of length 9 or longer can
generally be used to uniquely identify proteins.
