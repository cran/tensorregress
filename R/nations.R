#' nations data
#'
#' The array "R" is a 14 × 14 × 56 binary tensor consisting of 56 political relations of 14 countries between 1950 and 1965. The tensor entry indicates the presence or absence of a political action, such as “treaties”, “sends tourists to”, between the nations. Please set the diagonal elements Y(i,i,k) = 0 in the analysis.
#' The matrix "cov" is a 14 × 6 matrix describing a few important country attributes, e.g. whether a nation is actively involved in medicine NGO, law NGO, or belongs to a catholic nation, etc.
#'
#'@docType data
#'@usage data(nations)
#'
#'@format  A list. Includes a 14-14-56 binary array named "R" and a 14-6 matrix named "cov".
#'@keywords datasets
#'
#'
"nations"
