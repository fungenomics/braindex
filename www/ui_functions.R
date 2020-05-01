
library(shinycssloaders)

ws <- function(ui) withSpinner(ui, type = 3, size = 0.5,
                               color.background = "white", color = "#8896AE")

navigation <- function() {
	includeHTML('../www/layout/navigation.html')
}

beginPage <- function() {
  includeHTML('../www/layout/beginpage.html')
}

endPage <- function() {
  includeHTML('../www/layout/endpage.html')
}