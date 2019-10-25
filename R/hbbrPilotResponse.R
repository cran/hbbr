#' A list consisting of pilot data and associated discrete choice design information for the HBBR model framework.
#'
#' Data from 23 respondents each choosing preference from 18 choice cards.
#' Choice cards were randomly generated from 108 total choices. For details see article by
#' Mukhopadhyay, S., et al. "Hierarchical Bayesian Benefit-Risk Modeling
#' and Assessment Using Choice Based Conjoint." Statistics in Biopharmaceutical
#' Research 11.1 (2019): 52-60.
#'
#' @docType data
#'
#' @usage data(hbbrPilotResp)
#'
#' @format A list consisting of pilot response data from 23 experts and design information
#' to be used to fit the HBBR model
#' \describe{
#'    \item{brdta}{A data frame with 23x18 = 414 rows and 15 columns consists of responses from
#'    23 experts each providing tradeoff responses to 18 choice pairs. The 1st column consists of
#'    responders' id. The 2nd column contains binary responses (1 indicating 1st of the choice
#'    pair was selected, 0 indicating 2nd was selected). Remaining 13 columns
#'    contain the design matrix X taking values 0, 1, or -1; a value of 1 or -1 is used to
#'    indicate presence of an attribute level in the 1st choice or in the 2nd choice of the choice
#'    pair, respectively; a value of 0 is used to indicate absence of an attribute in the
#'    choice pair. See Details below for more about the discrete choice experiment that is coded as
#'    design matrix X.}
#'    \item{design}{A list of structure (b, r, bl, rl), where b and r indicate number of benefit
#'    and risk attributes, bl is a vector of integers of size b consisting number of levels within
#'    each benefit attribute; similarly rl is a vector of integers of size r consisting number
#'    of levels within each risk attribute.}
#' }
#' @details The discrete choice experiment (DCE) included 3 benefit
#'          attributes (b=3): overall survival (OS), objective response rate (ORR),
#'          fatigue reduction (FTG);
#'          and 2 risk attributes (r=2): febrile neutropenia (FebNEU) and severe pneumonia (SevPNA).
#'          There were 4 levels for each of the benefit attributes (ORR, OS, and FTG)
#'          (i.e. bl= rep(4,3)) and
#'          3 levels for each of the 2 risk attributes (FebNEU and SevPNA)
#'          (i.e. rl = rep(3,2)).
#'          The DCE produced b*r*(4 choose 2)*(3 choose 2) = 108 distinct non-dominant choice
#'          pairs each with one benefit and one
#'          risk attribute. Panels (questionnaires) were generated with 18 randomly selected
#'          choice pairs per panel from the set of 108 choice pairs.
#'          Since the part-worth of various levels within each attribute are to be measured
#'          relatively to the part-worth of the 1st level of the attribute, columns for
#'          the 1st level of the attributes are not required. Thus, we have
#'          sum(bl)-b + sum(br)-r = 13 columns are needed to
#'          obtain information on the X matrix which are stored as the last 13 columns of brdta.
#'
#' @keywords datasets
#'
#' @references Mukhopadhyay, S. et al. "Hierarchical Bayesian Benefit–Risk
#' Modeling and Assessment Using Choice Based Conjoint." Statistics in
#' Biopharmaceutical Research 11.1 (2019): 52-60.
#'
#' @examples
#' data(hbbrPilotResp)
"hbbrPilotResp"
