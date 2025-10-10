#' @title Shorten Each Word in a String
#'
#' @description
#' This function finds all words in an input string and replaces each word
#' with its first few letters. It can also remove specified substrings
#' before processing.
#'
#' @param input_str A character string to be processed.
#' @param num_letters (Optional) An integer specifying how many letters to
#'   keep from the beginning of each word. Must be 1 or greater. Default: 2.
#' @param stringsToRemove (Optional) A character vector of substrings to be
#'   completely removed from the input string before shortening. Default: `NULL`.
#'
#' @return A character string where each word has been replaced by its
#'   abbreviated version.
#'
#' @importFrom stringr str_remove_all str_replace_all
#'
#' @examples
#' # Example 1: Basic usage with default of 2 letters
#' original_text <- "This is a sample sentence for testing."
#' shorten_string(original_text)
#' # Expected output: "Th is a sa se fo te."
#'
#' # Example 2: Using 3 letters
#' shorten_string("Another example string", num_letters = 3)
#' # Expected output: "Ano exa str"
#'
#' # Example 3: Removing a specific string first
#' shorten_string("Remove this_word please", stringsToRemove = "this_word")
#' # Expected output: "Re  pl"
#' @export
shorten_string <- function(input_str, num_letters = 2, stringsToRemove = NULL) {
  # --- 1. Input Validation ---
  if (!is.numeric(num_letters) || num_letters < 1) {
    stop("`num_letters` must be a positive integer.")
  }
  num_letters <- as.integer(num_letters)

  # --- 2. Remove Specified Substrings ---
  if (!is.null(stringsToRemove)) {
    for (st in stringsToRemove) {
      input_str <- str_remove_all(input_str, st)
    }
  }

  # --- 3. Define the Replacement Function ---
  replacer_func <- function(word) {
    substr(word, 1, num_letters)
  }

  # --- 4. Perform the Replacement using stringr ---
  # The pattern `\\w+` matches one or more word characters.
  replaced_str <- str_replace_all(input_str, "\\w+", replacer_func)

  return(replaced_str)
}