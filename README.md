# Parallel-implementation-of-Sequence-Alignment-MPI-OpenMP-CUDA
Sequence Alignment – a way to estimate a similarity of two strings of letters - is an important field in bioinformatics .  Sequence is a string of capital letters including hyphen sign (-), for example

![image](https://user-images.githubusercontent.com/97550175/149371849-80886836-c78f-4246-8220-37092bf2251b.png)

Each letter in the sequence represents DNA, RNA, or protein. Identification of region of similarity of set of Sequences is extremely time consuming. 

Alignment Score Definition with pair-wise comparison

1.	Similarity of two sequences Seq1 and Seq2 of equal length is defined as follows:

•	Two Sequences are places one under another:

![image](https://user-images.githubusercontent.com/97550175/149371985-e5ed54d1-e961-41fe-b244-14135caca5d7.png)

•	Each letter from Seq1 is compared with the correspondent letter from Seq2. If these letters are identical the pair is marked with Star sign (*).

•	Otherwise the additional check is provided. The letters are checked if they both present at least in one of 9 groups called Conservative Groups:

![image](https://user-images.githubusercontent.com/97550175/149372051-6d535465-261f-41cc-8b04-2428551a094a.png)

In case that the pair is found in one of Conservative Group it is marked with Colon sign (:). 
For example, the pair (E, K) is marked with sign : because they both were found in group NEQK

•	If no Conservative Group is found, the pair is checked against 15 Semi-Conservative Groups

![image](https://user-images.githubusercontent.com/97550175/149372107-07dad46a-de8f-4f21-81ef-86535e4b192e.png)

If the pair do presents in one of Semi-Conservative Groups, it is marked with Point sign (.).
For example, the pair (K,S) is marked with sign . because they both were found in group STNK

•	If the letters in the pair are not equal, do not present both not in Conservative nor in Semi-Conservative groups – the pair is marked with Space sign (‘ ‘).

 
At the end of the check process the whole Sequence of Signs is obtained. This Sequence is used to estimate the similarity of two sequences – Seq1 and Seq2. For this project following formula is used to estimate the Alignment Score:
S = W1*NumberOfStars  -  W2*NumberOfColons – W3*NumberOfPoints  - W4*NumberOfSpaces
where Wi are the given weight coefficients.

2.	Similarity of two sequences Seq1 and Seq2 in case that Seq2 is shorter than Seq1, is defined as follows:

•	The Sequence Seq2 is places under the Sequence Seq1 with offset n from the start of the Sequence Seq1. The Sequence Seq2 do not allowed to pass behind the end of Seq1.
•	The letters from Seq1 that do not have a corresponding letter from Seq2 are ignored.
•	The Alignment Score is calculated according the pair-wise procedure described above.

For example, Sequence Seq2 is placed at different offsets under Seq1:
Score = -27
Offset n = 5
![image](https://user-images.githubusercontent.com/97550175/149372350-f64a353c-dd30-42e2-816c-a19c22556aae.png)

Score = 17
Offset n = 15
![image](https://user-images.githubusercontent.com/97550175/149372400-0fade5d5-2d49-472f-8034-c37e4705aa2e.png)

Mutant Sequence Definition

For a given Sequence S we define a Mutant Sequence MS(n) which is received by substitution of one or more letter by other letter defined by Substitution Rule. 
We will define the Substitution Rules as follows:
1.	The original letter is allowed to be substituted by another letter if there is no Conservative Group that contains both letters. For example, 
•	N is not allowed to be substituted by H because both letters present in Conservative Group NHQK
•	N may be substituted by W because there is now Conservative Group that contains both N and W
2.	It is not mandatory to substitute all instances of some letter by same substitution letter, for example the sequence  PSHLSPSQ has Mutant Sequence  PFHLSPLQ  

in this repository we will perform the following task - For two given sequences Seq1 and Seq2, find a mutant of Seq2 and its offset that produce a maximum / minimum Alignment Score.
