#Mads F. Schou
#Functions for data vizualization

#1. Pretty table for pdf knitting
#2. Pretty MCMCCglmmm table for pdf knitting

#######################################
###--- 1. Pretty table for pdf knitting
#######################################

prettytable = function(mydf, myfontsize = 9){
  kbl(mydf, booktabs = T, longtable = T) %>% 
    kable_styling(latex_options= c("HOLD_position","repeat_header"), position = "left",font_size = myfontsize) %>% 
    row_spec(0, bold = T) %>%
    column_spec(1, bold = T) %>%
    collapse_rows(columns = 1, latex_hline = "major", valign = "middle")
}

#######################################
###--- 2. Pretty MCMCCglmmm table for pdf knitting
#######################################

#See https://haozhu233.github.io/kableExtra/awesome_table_in_pdf.pdf
prettyMCtable = function(mydf, myfootnote = NULL, myfontsize = 7){
  #mydf_renam <- term.renaming(mydf,MyNewNames)
  mydf_renam <- mydf
  kbl(mydf_renam, booktabs = T) %>% 
    kable_styling(latex_options="HOLD_position", position = "left",font_size = myfontsize) %>% 
    row_spec(0, bold = T) %>%
    column_spec(1, bold = T) %>%
    column_spec(which(colnames(mydf_renam) == "pMCMC"),
                bold = as.numeric(as.character(mydf_renam[,"pMCMC"])) < "0.05" &
                  !is.na(as.numeric(as.character(mydf_renam[,"pMCMC"])))) %>%
    collapse_rows(columns = 1, latex_hline = "major", valign = "middle") %>%
    footnote(general = myfootnote, general_title = "")
}