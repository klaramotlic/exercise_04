#' Align two sequences globally using Hirschberg's algorithm
#' 
#' @param seq1 DNAString object representing NT or AA sequence to align
#' @param seq2 DNAString object representing NT or AA sequence to align
#' @param align A list of DNAString objects with alignment of input sequences
#' @param match An integer value of a score for matching bases
#' @param mismatch An integer value of score for mismatching bases
#' @param gap An integer value of penalty for gap insertion
hirschberg_template <- function(seq1, seq2, align, match, mismatch, gap){
    
    first_align_row <- align[[1]] # initialize the first row of alignment
    second_align_row <- align[[2]] # initialize the second row of alignment
  
  
    if (length(seq1) == 0)# length of seq1 is equal to zero
    {
        for (k in 1:length(seq2))# for each character in seq2
        {
            first_align_row <- append(first_align_row, '-')# add gap
            second_align_row <- append(second_align_row, seq2[k])# add character from seq2
        }
        align <- c(DNAStringSet(first_align_row), DNAStringSet(second_align_row))
        print(align)
    }
    else if (length(seq2) == 0)# length of seq2 is equal to zero
    {
        for (l in 1:length(seq1))# for each character in seq1
        {
            first_align_row <- append(first_align_row, seq1[l])# add character from seq1
            second_align_row <- append(second_align_row, '-')# add gap
        }
        align <- c(DNAStringSet(first_align_row), DNAStringSet(second_align_row))
        print(align)
    }
    else if (length(seq1) = 1 & length(seq2) = 1)# length of seq1 and seq2 is equal to 1
    {
        first_align_row <- append(first_align_row, seq1)# add character from seq1
        second_align_row <- append(second_align_row, seq2)# add character from seq2
        align <- c(DNAStringSet(first_align_row), DNAStringSet(second_align_row))
        print(align)
    }
    else
    {
        x_len <- length(seq1)# length of seq1
        x_mid <- floor(length(seq1)/2)# half of the length of seq1
        y_len <- length(seq2)# length of seq2
        
        left_score <- nw(seq1[1:x_mid], seq2)# NW score for the first half of seq1 and the whole seq2
        right_score <- nw(seq1[x_mid:end], seq2)# NW score for the second half of seq1 and the whole seq2 (both are reversed)
        y_mid <- max(left_score + rev(right_score))# index of division for seq2

        # The first half
        if (y_mid == 0)# index of division for seq2 is equal to 0
        {
            align <- hirshberg_template(seq1[1:x_mid], DNAString(), align, match, mismatch, gap)# call hirschberg function for the first half of seq1 and for an empty DNAString object
        }
        else
        {
            align <- hirshberg_template(seq1[1:x_mid], seq2[1:y_mid], align, match, mismatch, gap)# call hirschberg function for the first half of seq1 and for the first part of seq2
        }
        
        # The second half
        if ((x_mid + 1) > x_len) # seq1 cannot be further divided
        {
            align <- hirshberg_template(DNAString(), seq2[y_mid:end], align, match, mismatch, gap)# call hirschberg function for an empty DNAString object and the second half of seq2
        }
        else if ((y_mid + 1) > y_len) # seq2 cannot be further divided
        {
            align <- hirshberg_template(seq1[x_mid:end], DNAString(), align, match, mismatch, gap)# call hirschberg function for the second half of seq1 and for an empty DNAString object
        }
        else 
        {
            Align <- hirshberg_template(seq1[x_mid:end], seq2[y_mid:end], align, match, mismatch, gap)# call hirschberg function for the second half of seq1 and the second part of seq2
        }
    }
    return(align)
}


nw <- function(seq1, seq2){
  Score(0, 0) <- 0
  for (j in 1:length(seq2)){
    Score(0, j) <- Score(0, j - 1) + (seq2[j])
  }
  for (i in 1:length(seq1)){
    Score(1, 0) = Score(0, 0) + (seq1[i])
    for (j in 1:length(seq2)){
      scoreSub <- Score(0, j - 1) + (seq1[i], seq2[j])
      scoreDel <- Score(0, j) + (seq1[i])
      scoreIns <- Score(1, j - 1) + (seq2[j])
      Score(1, j) = max(scoreSub, scoreDel, scoreIns)
      end
    }
      
    }
  Score(0, :) = Score(1, :)
  for (j in 0:length(seq2)){
    LastLine[j] <- Score(1, j)
  }
  
  return LastLine
}