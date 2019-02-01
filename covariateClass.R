
setClass("Covariate",
         slots = c(
           Token_0 = "character",
           Token_1 = "character",
           Token_2 = "character",
           Token_3 = "character",
           Token_4 = "character",
           Token_5 = "character",
           interpolate = "logical"
         )
         )

#dictionaries
Token_0_dict <- c("chirps-v2-0", 
                     "LST_Day_v6",
                     "LST_Night_v6",
                     "TCB_v6.",
                     "EVI_v6.")

names(Token_0_dict) <- c("Rain",
                         "LST_Day",
                         "LST_Night",
                         "TCB",
                         "EVI") 

test <- new("Covariate", Token_0="1")

setGeneric("name", function(x) standardGeneric("name"))
setGeneric("name<-", function(x, value) standardGeneric("name<-"))

setMethod("name", "Covariate", function(x) names(Token_0_dict)[which(Token_0_dict==x@Token_0)])

setMethod("name<-", "Covariate", function(x, value) {
  x@Token_0 <- Token_0_dict[which(names(Token_0_dict) == value)]
  x
})


# Covariate <- function(name, 
#                       temporal,
#                       )

