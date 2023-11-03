package main

import (
	"fmt"
	"os"
	"github.com/TimothyStiles/poly/align"
	"github.com/TimothyStiles/poly/align/matrix"
	"github.com/TimothyStiles/poly/alphabet"
)

func main() {
	seq1 := os.Args[1]
	seq2 := os.Args[2]

	mat := [][]int{
                /*       A C G T U */
                /* A */ {1, -1, -1, -1, -1},
                /* C */ {-1, 1, -1, -1, -1},
                /* G */ {-1, -1, 1, -1, -1},
                /* T */ {-1, -1, -1, 1, -1},
                /* U */ {-1, -1, -1, -1, 1},
	}
	alphabet := alphabet.NewAlphabet([]string{"A", "C", "G", "T", "U"})
	
	subMatrix, err := matrix.NewSubstitutionMatrix(alphabet, alphabet, mat)
	if err != nil {
                fmt.Println("error: %s", err)
        }

	scoring, err := align.NewScoring(subMatrix, -1)
	if err != nil {
                fmt.Println("error: %s", err)
        }
	
	score, align1, align2, err := align.NeedlemanWunsch(seq1, seq2, scoring)
	if err != nil {
                fmt.Println("error: %s", err)
        }
	fmt.Printf("Score: %d\n", score)
	fmt.Println(align1)
	fmt.Println(align2)
}
