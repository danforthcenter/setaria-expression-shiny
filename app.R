list.of.packages <- c("shiny", "d3heatmap", "grid", 
                      "shinydashboard", "reshape", "reshape2", "stringr", 
                      "RColorBrewer","gridExtra","gplots","svglite",
                      "ggplot2", "gridExtra","shinyjs")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')

library(shiny)
library(ggplot2)
library(grid)
library(shinydashboard)
library(reshape)
library(reshape2)
library(stringr)
library(gridExtra)
library(gplots)
library(svglite)
library(d3heatmap)
library(shinyjs)

load_data <- function() {
  hide("loading_page")
  show("mainTabsetPanel")
}

########################################################################################
# Read Setaria Expression Data In (See methods for details on generation of datatables)
# See setaria-rnaseq-prep.R script for how data was cleaned up / normalized
########################################################################################

experiment<-readRDS('data/experiment.information.rds')
setaria.data<-readRDS('data/setaria-all-data.rds')

########################################################################################
# Setaria Shiny Application -ui.R
########################################################################################
ui <- dashboardPage(

  dashboardHeader(title="Setaria Expression Explorer"),

  dashboardSidebar(disable = T),
  
  dashboardBody(
    useShinyjs(),
    div(
      id = "loading_page",
      h1("Loading...")
    ),
    hidden(
      div(id="mainTabsetPanel",
        tabsetPanel(
          
          ########################################################################################
          tabPanel(title="Welcome",
               fluidRow(
                 tags$head(includeScript("google-analytics.js")),
                 box(title="Setaria Expression Explorer",width=12, solidHeader = T,status = 'primary',
                     p("This tool is brought to you by the",a(href="http://www.gehan-lab.org/",target='_blank',"Gehan Lab"),
                       "at the Donald Danforth Plant Science Center. For the code used to generate this app,
                       please visit ", a(href="https://github.com/danforthcenter/setaria-expression-shiny",target='_blank',"our Github"),
                       ". For more information on the experiments please refer to the primary sources: For the leaf data see ", 
                       a(href="https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1080-3",target='_blank',"Studer et al. 2016"),
                       " and for the infloresence data see ", 
                       a(href="https://www.frontiersin.org/articles/10.3389/fpls.2018.01309/full#h3",target='_blank',"Zhu et al. 2018"))
                     ),
                 
                 box(title="Using this Tool",width=12, solidHeader = T,status = 'primary',
                      p("In the SAMPLE INFO tab see the available datasets and conditions"),
                      hr(),
                      p("In the SEARCH AND BROWSE DATA tab above, you can search for genes of interest using  by either
                      the search bar or by uploading a .txt file with gene ids. Alternatively, you can use the data filters
                      to browse the data.Once you have loaded or searched for your selections the data can be viewed and dowloaded. 
                      A plot of the data can be seen on the Plot Data tab once the plot button is hit."),
                      hr(),
                      p("In the PLOT DATA tab above, data selected from Select Data or Browse Data tabs is plotted")
                    )
                   )
                  ),
          ########################################################################################
          
          tabPanel(title="Sample Info",
                   fluidRow(
                     box(title="Sample Info",width=12, solidHeader = T,status = 'danger',
                         p("For more information on the experiments please refer to the primary data sources: 
                           For the leaf data see ", 
                           a(href="https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1080-3",target='_blank',"Studer et al. 2016")," and for the infloresence data see ", 
                           a(href="https://www.frontiersin.org/articles/10.3389/fpls.2018.01309/full#h3",target='_blank',"Zhu et al. 2018")),
                         div(style = 'overflow-x: scroll', dataTableOutput("conditions.table"))
                     )
                   )
          ),
          ########################################################################################
          
          tabPanel(title="Search and Browse Data",
                   fluidRow(
                     box(title="Search Data with GENEID or GO",width=12, solidHeader = T,status = 'success',collapsible = TRUE, collapsed = TRUE,
                         h4("Search using small sets of GENEIDs"),
                         p("GENEIDs, Orthologs, or GO separated by a comma are allowed"),
                         textInput('gene_search_text', label="example: Sevir.1G000100", value = ""),
                         h4("Search using small sets of GO TERMS (Pfam,Panther,KOG,KEGG,KO,GO)"),
                         textInput('go_search_text', label="example:GO:GO:0008270", value = ""),
                         h4("Search using small sets of ORTHOLOG GENEIDs (S. italica, Arabidiopsis, Rice, Brachypodium, Maize, Sorghum)"),
                         textInput('orth_search_text', label="example:AT5G19850.1,AT5G26667.3", value = ""),
                         p("refresh page to clear search")
                     )),
                   
                   fluidRow(
                     box(title="Search Data with File",width=12, solidHeader = T,status = 'success',collapsible = TRUE,collapsed=TRUE,
                         h4("Upload a file of GENEIDs, the GENEIDs should not be quoted, one line per geneid"),
                         fileInput('file.geneid', 'Choose file to upload',
                                   accept = c('text/csv','text/comma-separated-values','text/tab-separated-values','text/plain', '.csv','.tsv')),
                         checkboxInput('header', 'Header', FALSE),
                         h4("Upload a file of ORTHOLOG GENEIDs, the ORTHOLOG GENEIDs should not be quoted, one line per geneid"),
                         fileInput('file.ortholog', 'Choose file to upload',
                                   accept = c('text/csv','text/comma-separated-values','text/tab-separated-values','text/plain', '.csv','.tsv')),
                         checkboxInput('header1', 'Header', FALSE),
                         p("refresh page to clear search")
                     )),
                   fluidRow(
                     box(title="Genes, Orthologs, or GO Selected with Search",width=12, solidHeader = T,status = 'success',collapsible = TRUE,collapsed =TRUE,
                         verbatimTextOutput("selected.genes"),
                         verbatimTextOutput("selected.orthologs"),
                         verbatimTextOutput("selected.go")
                         )),
                   fluidRow(
                     box(title="Filter Data by Minimum TPM",width=12, solidHeader = T,status = 'success', collapsible = TRUE,collapsed =TRUE,
                         fluidRow(
                                column(6,
                                       radioButtons('cuttype', 'Cut-Off Type',c('Meets Cut-Off In All Leaf Samples'='allleaf', 
                                                                                'Meets Cut-Off In All Inflorescence Samples'='allinflor',
                                                                                'Meets Cut-Off All Leaf and Inflorescence Samples'='alls',
                                                                                'Meets Cut-Off In At Least One Leaf Sample'='oneleaf',
                                                                                'Meets Cut-Off In At Least One Inflorescence Sample'='oneinflor',
                                                                                'Meets Cut-Off In At Least One Leaf or Inflorescence Samples'='onesamp'
                                                                                ),'alls')),
                                
                                column(6,
                                       numericInput("tpm.cutoff",
                                                   "Minimum TPM Cutoff (0 to 66468):",
                                                   0, min = 0, max = 66468),
                                                   verbatimTextOutput("value"))

                      ))),
                   fluidRow(
                     box(title="Filter Data by SD (percent of mean)",width=12, solidHeader = T,status = 'success', collapsible = TRUE,collapsed =TRUE,
                         fluidRow(
                           column(6,
                                  radioButtons('cuttype1', 'Cut-Off Type',c('Meets Cut-Off In All Leaf Samples'='allleaf', 
                                                                           'Meets Cut-Off In All Inflorescence Samples'='allinflor',
                                                                           'Meets Cut-Off All Leaf and Inflorescence Samples'='alls',
                                                                           'Meets Cut-Off In At Least One Leaf Sample'='oneleaf',
                                                                           'Meets Cut-Off In At Least One Inflorescence Sample'='oneinflor',
                                                                           'Meets Cut-Off In At Least One Leaf or Inflorescence Samples'='onesamp'
                                  ),'alls')),
                           
                           column(6,
                                  sliderInput("sd.cutoff",
                                              "Percent Standard Deviation Cut-Off (%):",
                                              min = 0, max = 100,
                                              value = 100))
                           
                         ))),
                   fluidRow(
                     box(title="Data",width=12, solidHeader = T,status = 'success', collapsible = TRUE,collapsed =FALSE,
                           div(style = 'overflow-x: scroll', dataTableOutput("setaria.data"))
                     )),
                     fluidRow(
                       box(title="Download Selected Data",width=12, solidHeader = T,status = 'success',
                             column(4,
                                    downloadButton('download.selected', "Download Selected Data"))
                           ))
                     ),
          ########################################################################################
          tabPanel(title="Plot Data",
                   fluidRow(
                     box(title="Plot Data",width=12, solidHeader = T,status = 'info',
                        h4(textOutput("numbergenes")),
                        h4("We restrict plotting to 100 genes. 
                            Please download data from 'Search and Browse Data' tab or run a local installation of", 
                            a(href="https://github.com/maliagehan/diel-explorer/",target='_blank',"Setaria Expression Explorer"), 
                            "if you want to graph more genes."))),
                   fluidRow(
                     box(title="Plot Line Graph",width=12, solidHeader = T,status = 'info',collapsible = TRUE,collapsed =TRUE,
                         actionButton("plot.line", label="Plot Selected Data as Line Graph" ),
                         h4("Leaf Gradient Data"),
                         plotOutput("leaf.line.plot"),
                         downloadButton('download.leaf.plot',"Download Leaf Plot"),
                         h4("Infloresence Development"),
                         plotOutput("inflor.line.plot"),
                         downloadButton('download.inflor.plot',"Download Infloresence Plot"))),
                     fluidRow(
                       box(title="Plot Heatmap",width=12, solidHeader = T,status = 'info',collapsible = TRUE,collapsed =TRUE,
                         actionButton("plot.heat", label="Plot Selected Data as Heatmap" ),
                         radioButtons('rowcol1', 'Scale Heatmap',c(Row='TRUE',Column='FALSE'),'TRUE'),
                         h4("Leaf Gradient Data"),
                         d3heatmapOutput("leaf.heatmap"),
                         downloadButton('download.leaf.heat',"Download Leaf Heatmap"),
                         h4("Infloresence Development"),
                         d3heatmapOutput("inflor.heatmap"),
                         downloadButton('download.inflor.heat',"Download Infloresence Heatmap")))
                 ),
          ########################################################################################
          tabPanel(title="Contact Us",
                   fluidRow(
                     box(title="Contact US",width=12, solidHeader = T,status = 'danger',
                         p("For questions or if you are interested in adding your own data contact ",
                           a(href="https://github.com/danforthcenter/setaria-expression-shiny/issues",target='_blank',"Malia Gehan")),
                         p("For more information on the Gehan lab visit our ",
                           a(href="http://www.gehan-lab.org/",target='_blank',"website"))
                     )
                   )
                )
        ))
    )
  )
)
    
########################################################################################
# Setaria Shiny Application -server.R
########################################################################################
server<-function(input,output,session){
  
  #output for sample info tab ########################################################
    
  output$conditions.table<-renderDataTable(experiment, options=list(paging=FALSE,searching=FALSE))
  
  #output for select data tab #########################################################
  
  #output for browse data tab #########################################################
  
  #get the search terms
  searchgenes<-reactive({
    genes<-gsub(" ","",input$gene_search_text)
    genes1<-data.frame(strsplit(genes,","))
    colnames(genes1)<-c("GENEID")
    genes1
    })
  
  #get the file contents
  searchfile<-reactive({
    geneids<-input$file.geneid
    
    if(is.null(geneids))
      return(NULL)
    
    genes1<-read.csv(geneids$datapath,header=input$header, strip.white = TRUE)
    colnames(genes1)<-c("GeneID")
    genes1
  })
  
  #join the gene search and the file contents
  joinedsearch<-reactive({
    rbind(searchgenes(),searchfile())
  })
  
  #output the list of genes so the user can see it
  output$selected.genes<-renderPrint({
    joinedsearch()
  })
  
  #get the search ortholog terms
  searchorthologs<-reactive({
    orth<-gsub(" ","",input$orth_search_text)
    orth1<-data.frame(strsplit(orth,","))
    colnames(orth1)<-c("ORTHOLOGS")
    orth1
  })
  
  #get the ortholog file terms
  searchfileorth<-reactive({
    orth<-input$file.ortholog
    
    if(is.null(orth))
      return(NULL)
    
    orth1<-read.csv(orth$datapath,header=input$header1,strip.white =TRUE)
    colnames(orth1)<-c("ORTHOLOGS")
    orth1
  })
  
  #join the ortholog search and the ortholog file contents
  joinedsearchorth<-reactive({
    rbind(searchorthologs(),searchfileorth())
  })
  
  #output the list of genes so the user can see it
  output$selected.orthologs<-renderPrint({
    joinedsearchorth()
  })
  
  #get list of go terms
  searchgo<-reactive({
    go<-gsub(" " ,"",input$go_search_text)
    go1<-data.frame(strsplit(go,","))
    colnames(go1)<-c("GO-TERM")
    go1
  })
  
  #output the list of go so the user can see it
  output$selected.go<-renderPrint({
    searchgo()
  })
  
  #filter data
  setaria.input<-reactive({
    
    setaria.table<-setaria.data
    
    if(nrow(joinedsearch())!=0 | nrow(searchgo())!=0 | nrow(joinedsearchorth())!=0){
       colnum<-ncol(setaria.table)
       setaria.subset <- data.frame(matrix(ncol=colnum, nrow = 0))
       colnames(setaria.subset) <- paste0(c(colnames(setaria.table)))
    
       if(nrow(joinedsearch())!=0){
         for(x in 1:nrow(joinedsearch())){
           search<-as.character(joinedsearch()[x,1])
           row<-setaria.table[(grep(search,setaria.table$GeneID)),]
           setaria.subset<-rbind(row,setaria.subset)
         }
       }
       
       if(nrow(joinedsearchorth())!=0){
         for(x in 1:nrow(joinedsearchorth())){
           search<-as.character(joinedsearchorth()[x,1])
           row1<-setaria.table[(grep(search,setaria.table$Best.hit.arabi.name)),]
           row2<-setaria.table[(grep(search,setaria.table$Best.hit.rice.name)),]
           row3<-setaria.table[(grep(search,setaria.table$Maizev4.gene)),]
           row4<-setaria.table[(grep(search,setaria.table$Zmays_284.gene)),]
           row5<-setaria.table[(grep(search,setaria.table$Bradi.gene)),]
           row6<-setaria.table[(grep(search,setaria.table$Sobic.gene)),]
           row7<-setaria.table[(grep(search,setaria.table$Seita.gene)),]
           setaria.subset<-rbind(row1,setaria.subset)
           setaria.subset<-rbind(row2,setaria.subset)
           setaria.subset<-rbind(row3,setaria.subset)
           setaria.subset<-rbind(row4,setaria.subset)
           setaria.subset<-rbind(row5,setaria.subset)
           setaria.subset<-rbind(row6,setaria.subset)
           setaria.subset<-rbind(row7,setaria.subset)
         }
       }
    
       if(nrow(searchgo())!=0){
         for(x in 1:nrow(searchgo())){
           search<-as.character(searchgo()[x,1])
           row1<-setaria.table[(grep(search,setaria.table$Pfam)),]
           row2<-setaria.table[(grep(search,setaria.table$Panther)),]
           row3<-setaria.table[(grep(search,setaria.table$KOG)),]
           row4<-setaria.table[(grep(search,setaria.table$KEGG.ec)),]
           row5<-setaria.table[(grep(search,setaria.table$KO)),]
           row6<-setaria.table[(grep(search,setaria.table$GO)),]
           setaria.subset<-rbind(row1,setaria.subset)
           setaria.subset<-rbind(row2,setaria.subset)
           setaria.subset<-rbind(row3,setaria.subset)
           setaria.subset<-rbind(row4,setaria.subset)
           setaria.subset<-rbind(row5,setaria.subset)
           setaria.subset<-rbind(row6,setaria.subset)
           
         }
       }
       
       setaria.table<-setaria.subset
    }
    
    if(input$cuttype=="allleaf"){
      setaria.table<-setaria.table
      if (input$tpm.cutoff!=0){
        setaria.table5<-setaria.table[setaria.table$means1>=as.numeric(input$tpm.cutoff),]
        setaria.table4<-setaria.table5[setaria.table5$means2>=as.numeric(input$tpm.cutoff),]
        setaria.table3<-setaria.table4[setaria.table4$means3>=as.numeric(input$tpm.cutoff),]
        setaria.table2<-setaria.table3[setaria.table3$means4>=as.numeric(input$tpm.cutoff),]
        setaria.table1<-setaria.table2[rowSums(is.na(setaria.table2)) != ncol(setaria.table2), ]
        setaria.table<-setaria.table1
      }
    }
    
    if(input$cuttype=="allinflor"){
      setaria.table<-setaria.table
      if (input$tpm.cutoff!=0){
        setaria.table7<-setaria.table[setaria.table$X10DAS_mean>=as.numeric(input$tpm.cutoff),]
        setaria.table6<-setaria.table7[setaria.table7$X12DAS_mean>=as.numeric(input$tpm.cutoff),]
        setaria.table5<-setaria.table6[setaria.table6$X14DAS_mean>=as.numeric(input$tpm.cutoff),]
        setaria.table4<-setaria.table5[setaria.table5$X15DAS_mean>=as.numeric(input$tpm.cutoff),]
        setaria.table3<-setaria.table4[setaria.table4$X16DAS_mean>=as.numeric(input$tpm.cutoff),]
        setaria.table2<-setaria.table3[setaria.table3$X18DAS_mean>=as.numeric(input$tpm.cutoff),]
        setaria.table1<-setaria.table2[rowSums(is.na(setaria.table2)) != ncol(setaria.table2), ]
        setaria.table<-setaria.table1
      }
    }

    if(input$cuttype=="alls"){
      if (input$tpm.cutoff!=0){
        setaria.table11<-setaria.table[setaria.table$means1>=as.numeric(input$tpm.cutoff),]
        setaria.table10<-setaria.table11[setaria.table11$means2>=as.numeric(input$tpm.cutoff),]
        setaria.table9<-setaria.table10[setaria.table10$means3>=as.numeric(input$tpm.cutoff),]
        setaria.table8<-setaria.table9[setaria.table9$means4>=as.numeric(input$tpm.cutoff),]
        setaria.table7<-setaria.table8[setaria.table8$X10DAS_mean>=as.numeric(input$tpm.cutoff),]
        setaria.table6<-setaria.table7[setaria.table7$X12DAS_mean>=as.numeric(input$tpm.cutoff),]
        setaria.table5<-setaria.table6[setaria.table6$X14DAS_mean>=as.numeric(input$tpm.cutoff),]
        setaria.table4<-setaria.table5[setaria.table5$X15DAS_mean>=as.numeric(input$tpm.cutoff),]
        setaria.table3<-setaria.table4[setaria.table4$X16DAS_mean>=as.numeric(input$tpm.cutoff),]
        setaria.table2<-setaria.table3[setaria.table3$X18DAS_mean>=as.numeric(input$tpm.cutoff),]
        setaria.table1<-setaria.table2[rowSums(is.na(setaria.table2)) != ncol(setaria.table2), ]
        setaria.table<-setaria.table1
      }
    }

    if(input$cuttype=="oneleaf"){
      if (input$tpm.cutoff!=0){
        setaria.table2<-setaria.table[(setaria.table$means1>=as.numeric(input$tpm.cutoff) |
                                       setaria.table$means2>=as.numeric(input$tpm.cutoff) |
                                       setaria.table$means3>=as.numeric(input$tpm.cutoff)|
                                       setaria.table$means4>=as.numeric(input$tpm.cutoff)),]
        setaria.table1<-setaria.table2[rowSums(is.na(setaria.table2)) != ncol(setaria.table2), ]
        setaria.table<-setaria.table1
      }
    }

    if(input$cuttype=="oneinflor"){
      if (input$tpm.cutoff!=0){
        setaria.table2<-setaria.table[(setaria.table$X10DAS_mean>=as.numeric(input$tpm.cutoff) |
                                         setaria.table$X12DAS_mean>=as.numeric(input$tpm.cutoff) |
                                         setaria.table$X14DAS_mean>=as.numeric(input$tpm.cutoff)|
                                         setaria.table$X15DAS_mean>=as.numeric(input$tpm.cutoff)|
                                         setaria.table$X16DAS_mean>=as.numeric(input$tpm.cutoff)|
                                         setaria.table$X18DAS_mean>=as.numeric(input$tpm.cutoff)),]
        setaria.table1<-setaria.table2[rowSums(is.na(setaria.table2)) != ncol(setaria.table2), ]
        setaria.table<-setaria.table1
      }
    }

    if(input$cuttype=="onesamp"){
      if (input$tpm.cutoff!=0){
        setaria.table6<-setaria.table[(setaria.table$means1>=as.numeric(input$tpm.cutoff) |
                                         setaria.table$means2>=as.numeric(input$tpm.cutoff) |
                                         setaria.table$means3>=as.numeric(input$tpm.cutoff)|
                                         setaria.table$means4>=as.numeric(input$tpm.cutoff) |
                                         setaria.table$X10DAS_mean>=as.numeric(input$tpm.cutoff) |
                                         setaria.table$X12DAS_mean>=as.numeric(input$tpm.cutoff) |
                                         setaria.table$X14DAS_mean>=as.numeric(input$tpm.cutoff)|
                                         setaria.table$X15DAS_mean>=as.numeric(input$tpm.cutoff)|
                                         setaria.table$X16DAS_mean>=as.numeric(input$tpm.cutoff)|
                                         setaria.table$X18DAS_mean>=as.numeric(input$tpm.cutoff)),]
        setaria.table5<-setaria.table6[rowSums(is.na(setaria.table6)) != ncol(setaria.table6), ]
        setaria.table<-setaria.table5
      }
    }

    if(input$cuttype1=="allleaf"){
      setaria.table<-setaria.table
      if (input$sd.cutoff!=100){
        setaria.table5<-setaria.table[setaria.table$percents1<=as.numeric(input$sd.cutoff),]
        setaria.table4<-setaria.table5[setaria.table5$percents2<=as.numeric(input$sd.cutoff),]
        setaria.table3<-setaria.table4[setaria.table4$percents3<=as.numeric(input$sd.cutoff),]
        setaria.table2<-setaria.table3[setaria.table3$percents4<=as.numeric(input$sd.cutoff),]
        setaria.table1<-setaria.table2[rowSums(is.na(setaria.table2)) != ncol(setaria.table2), ]
        setaria.table<-setaria.table1
      }
    }
    
    if(input$cuttype1=="allinflor"){
      setaria.table<-setaria.table
      if (input$sd.cutoff!=100){
        setaria.table7<-setaria.table[setaria.table$percent10d<=as.numeric(input$sd.cutoff),]
        setaria.table6<-setaria.table7[setaria.table7$percent12d<=as.numeric(input$sd.cutoff),]
        setaria.table5<-setaria.table6[setaria.table6$percent14d<=as.numeric(input$sd.cutoff),]
        setaria.table4<-setaria.table5[setaria.table5$percent15d<=as.numeric(input$sd.cutoff),]
        setaria.table3<-setaria.table4[setaria.table4$percent16d<=as.numeric(input$sd.cutoff),]
        setaria.table2<-setaria.table3[setaria.table3$percent18d<=as.numeric(input$sd.cutoff),]
        setaria.table1<-setaria.table2[rowSums(is.na(setaria.table2)) != ncol(setaria.table2), ]
        setaria.table<-setaria.table1
      }
    }
    
    if(input$cuttype1=="alls"){
      if (input$sd.cutoff!=100){
        setaria.table11<-setaria.table[setaria.table$percents1<=as.numeric(input$sd.cutoff),]
        setaria.table10<-setaria.table11[setaria.table11$percents2<=as.numeric(input$sd.cutoff),]
        setaria.table9<-setaria.table10[setaria.table10$percents3<=as.numeric(input$sd.cutoff),]
        setaria.table8<-setaria.table9[setaria.table9$percents4<=as.numeric(input$sd.cutoff),]
        setaria.table7<-setaria.table8[setaria.table8$percent10d<=as.numeric(input$sd.cutoff),]
        setaria.table6<-setaria.table7[setaria.table7$percent12d<=as.numeric(input$sd.cutoff),]
        setaria.table5<-setaria.table6[setaria.table6$percent14d<=as.numeric(input$sd.cutoff),]
        setaria.table4<-setaria.table5[setaria.table5$percent15d<=as.numeric(input$sd.cutoff),]
        setaria.table3<-setaria.table4[setaria.table4$percent16d<=as.numeric(input$sd.cutoff),]
        setaria.table2<-setaria.table3[setaria.table3$percent18d<=as.numeric(input$sd.cutoff),]
        setaria.table1<-setaria.table2[rowSums(is.na(setaria.table2)) != ncol(setaria.table2), ]
        setaria.table<-setaria.table1
      }
    }
    
    if(input$cuttype1=="oneleaf"){
      if (input$sd.cutoff!=100){
        setaria.table2<-setaria.table[(setaria.table$percents1<=as.numeric(input$sd.cutoff) |
                                         setaria.table$percents2<=as.numeric(input$sd.cutoff) |
                                         setaria.table$percents3<=as.numeric(input$sd.cutoff)|
                                         setaria.table$percents4<=as.numeric(input$sd.cutoff)),]
        setaria.table1<-setaria.table2[rowSums(is.na(setaria.table2)) != ncol(setaria.table2), ]
        setaria.table<-setaria.table1
      }
    }
    
    if(input$cuttype1=="oneinflor"){
      if (input$sd.cutoff!=100){
        setaria.table2<-setaria.table[(setaria.table$percent10d<=as.numeric(input$sd.cutoff) |
                                         setaria.table$percent12d<=as.numeric(input$sd.cutoff) |
                                         setaria.table$percent14d<=as.numeric(input$sd.cutoff)|
                                         setaria.table$percent15d<=as.numeric(input$sd.cutoff)|
                                         setaria.table$percent16d<=as.numeric(input$sd.cutoff)|
                                         setaria.table$percent18d<=as.numeric(input$sd.cutoff)),]
        setaria.table1<-setaria.table2[rowSums(is.na(setaria.table2)) != ncol(setaria.table2), ]
        setaria.table<-setaria.table1
      }
    }
    
    if(input$cuttype1=="onesamp"){
      if (input$sd.cutoff!=100){
        setaria.table6<-setaria.table[(setaria.table$percents1<=as.numeric(input$sd.cutoff) |
                                         setaria.table$percents2<=as.numeric(input$sd.cutoff) |
                                         setaria.table$percents3<=as.numeric(input$sd.cutoff)|
                                         setaria.table$percents4<=as.numeric(input$sd.cutoff) |
                                         setaria.table$percent10d<=as.numeric(input$sd.cutoff) |
                                         setaria.table$percent12d<=as.numeric(input$sd.cutoff) |
                                         setaria.table$percent14d<=as.numeric(input$sd.cutoff)|
                                         setaria.table$percent15d<=as.numeric(input$sd.cutoff)|
                                         setaria.table$percent16d<=as.numeric(input$sd.cutoff)|
                                         setaria.table$percent18d<=as.numeric(input$sd.cutoff)),]
        setaria.table5<-setaria.table6[rowSums(is.na(setaria.table6)) != ncol(setaria.table6), ]
        setaria.table<-setaria.table5
      }
    }
    setaria.table
  })
    
  output$setaria.data<-renderDataTable(({
    setaria.input()}),options=list(searching=FALSE, na="NA"))
  
  output$download.selected <- downloadHandler(
    filename = function() { paste('setaria_data_',strftime(Sys.time(),"%m-%d-%y_%H%M%S"),'.csv', sep='') },
    content = function(file) {
      write.csv(setaria.input(), file)
      }
    )

  #output for plot data tab ###########################################################  
  output$numbergenes<-reactive({
    paste(nrow(setaria.input()), "Genes are currently selected")
  })
  

  plot.input1<-reactive({

    setaria.graph<-setaria.input()

    if(nrow(setaria.graph)>100){
      setaria.plot1=setaria.graph[1:100,]
    }else{setaria.plot1=setaria.graph}

    setaria.plot<-subset(setaria.plot1,select=c(means1, means2,means3,means4))
    rownames(setaria.plot)<-setaria.plot1$GeneID
    return(setaria.plot)
  })
  
  plot.input3<-reactive({
    
    setaria.graph<-setaria.input()
    
    if(nrow(setaria.graph)>100){
      setaria.plot1=setaria.graph[1:100,]
    }else{setaria.plot1=setaria.graph}
    
    setaria.plot<-subset(setaria.plot1,select=c(X10DAS_mean,X12DAS_mean,X14DAS_mean,X15DAS_mean,X16DAS_mean,X18DAS_mean))
    rownames(setaria.plot)<-setaria.plot1$GeneID
    return(setaria.plot)
    })
  
  
  plot.input5<-reactive({
    setaria.graph<-setaria.input()
    
    if(nrow(setaria.graph)>100){
      setaria.plot1=setaria.graph[1:100,]
    }else{setaria.plot1=setaria.graph}
    
    setaria.plot<-subset(setaria.plot1,select=c(GeneID, means1, means2,means3,means4))
    setaria.melt<-melt(setaria.plot,by=GeneID)
    colnames(setaria.melt)<-c("GeneID","leaf.segment","mean")
    
    setaria.sd<-subset(setaria.plot1,select=c(GeneID, sd1,sd2,sd3,sd4))
    colnames(setaria.sd)<-c("GeneID","means1","means2","means3","means4")                  
    setaria.sd.melt<-melt(setaria.sd,by=GeneID)
    colnames(setaria.sd.melt)<-c("GeneID","leaf.segment","sd")
    
    setaria.plot.sd<-merge(setaria.melt,setaria.sd.melt,by=c("GeneID","leaf.segment"))

    ggplot(setaria.plot.sd, aes(x=leaf.segment, y=mean, group=factor(GeneID), colour=factor(GeneID))) +
      geom_line()+
      geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) +
      geom_point()+
      theme_bw()+
      labs(title="Leaf Segment Expression", y="TPM", x="Leaf Segment")+
      theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right")
  })
  
  plot.input6<-reactive({
    setaria.graph<-setaria.input()
    
    if(nrow(setaria.graph)>100){
      setaria.plot1=setaria.graph[1:100,]
    }else{setaria.plot1=setaria.graph}
    
    setaria.plot<-subset(setaria.plot1,select=c(GeneID, X10DAS_mean,X12DAS_mean,X14DAS_mean,X15DAS_mean,X16DAS_mean,X18DAS_mean))
    setaria.melt<-melt(setaria.plot,by=GeneID)
    colnames(setaria.melt)<-c("GeneID","leaf.segment","mean")
    
    setaria.sd<-subset(setaria.plot1,select=c(GeneID, X10DAS_SD,X12DAS_SD,X14DAS_SD,X15DAS_SD,X16DAS_SD,X18DAS_SD))
    colnames(setaria.sd)<-c("GeneID","X10DAS_mean","X12DAS_mean","X14DAS_mean","X15DAS_mean","X16DAS_mean","X18DAS_mean")                  
    setaria.sd.melt<-melt(setaria.sd,by=GeneID)
    colnames(setaria.sd.melt)<-c("GeneID","leaf.segment","sd")
    
    setaria.plot.sd<-merge(setaria.melt,setaria.sd.melt,by=c("GeneID","leaf.segment"))

    ggplot(setaria.plot.sd, aes(x=leaf.segment, y=mean, group=factor(GeneID), colour=factor(GeneID))) +
      geom_line()+
      geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) +
      geom_point()+
      theme_bw()+
      labs(title="Infloresence Development", y="TPM", x="Developmental Stage")+
      theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right")
  })
  
  plot.input2<-reactive({
    setaria.matrix=plot.input3()
    output$inflor.heatmap<-renderD3heatmap({
        if(input$rowcol1=='FALSE'){
          color.palette = colorRampPalette(c("lightyellow","lightblue","blue"),space="rgb")
          d3heatmap(setaria.matrix, scale="column",dendrogram = "none",margins=c(40, 200), color=color.palette(256))
        }else{
          color.palette = colorRampPalette(c("lightyellow","lightblue","blue"),space="rgb")
          d3heatmap(setaria.matrix, scale="row",dendrogram = "none",margins=c(40, 200),color=color.palette(256))}
      })
  })
  
  
  plot.input4<-reactive({
    setaria.matrix=plot.input1()
    output$leaf.heatmap<-renderD3heatmap({
      if(input$rowcol1=='FALSE'){
        color.palette = colorRampPalette(c("lightyellow","lightblue","blue"),space="rgb")
        d3heatmap(setaria.matrix, scale="column",dendrogram = "none",margins=c(40, 200), color=color.palette(256))
      }else{
        color.palette = colorRampPalette(c("lightyellow","lightblue","blue"),space="rgb")
        d3heatmap(setaria.matrix, scale="row",dendrogram = "none",margins=c(40, 200),color=color.palette(256))
      }})
  })

  observeEvent(input$plot.line,{
    output$leaf.line.plot<-renderPlot({
      withProgress(message = 'Making plot', value =NULL, {
        incProgress()
        grid.arrange(plot.input5(),ncol=1)
        })
    })
    output$inflor.line.plot<-renderPlot({
      withProgress(message = 'Making plot', value =NULL, {
        incProgress()
        grid.arrange(plot.input6(),ncol=1)
      })
    })
    })
  
  observeEvent(input$plot.line,{
    withProgress(message='Making plot', value=NULL,{
      incProgress()
      plot.input5()
    })})
  
  observeEvent(input$plot.heat,{
    withProgress(message='Making Heatmap', value=NULL,{
    incProgress()
    plot.input4()
    plot.input2()
    })})
  
  output$download.leaf.plot <- downloadHandler(
    filename = function() { paste('setaria_leaf_plot_',strftime(Sys.time(),"%m-%d-%y_%H%M%S"),".svg", sep='') },
    content=function(file){
      ggsave(file, plot = plot.input1(), device = "svg",height =8, width=10)
    }, contentType='image/svg')
  
  output$download.inflor.plot <- downloadHandler(
    filename = function() { paste('setaria_inflor_plot_',strftime(Sys.time(),"%m-%d-%y_%H%M%S"),".svg", sep='') },
    content=function(file){
      ggsave(file, plot = plot.input3(), device = "svg",height =8, width=10)
    }, contentType='image/svg')
  
  output$download.leaf.heat <- downloadHandler(
    filename = function() { paste('setaria_leaf_heatmap_',strftime(Sys.time(),"%m-%d-%y_%H%M%S"),".pdf", sep='') },
    content=function(file){
      color.palette = colorRampPalette(c("lightyellow","lightblue","blue"),space="rgb")
      if(input$rowcol1=='FALSE'){scale1="column"}else{scale1="none"}
      pdf(file, width=8, height=7, useDingbats =FALSE)
      heatmap.2(as.matrix(plot.input1()),
                Rowv=FALSE,
                Colv=FALSE,
                dendrogram='none',
                scale=scale1,
                col=color.palette(256),
                trace='none',
                margins=c(8,20),
                symbreaks=FALSE,
                symm=FALSE,
                cex.main=0.75,
                density.info="none"
      )
      dev.off()
    },contentType='image/pdf')
  
  output$download.inflor.heat <- downloadHandler(
    filename = function() { paste('setaria_inflor_heatmap_',strftime(Sys.time(),"%m-%d-%y_%H%M%S"),".pdf", sep='') },
    content=function(file){
      color.palette = colorRampPalette(c("lightyellow","lightblue","blue"),space="rgb")
      if(input$rowcol1=='FALSE'){scale1="column"}else{scale1="none"}
      pdf(file, width=8, height=7, useDingbats =FALSE)
      heatmap.2(as.matrix(plot.input3()),
                Rowv=FALSE,
                Colv=FALSE,
                dendrogram='none',
                scale=scale1,
                col=color.palette(256),
                trace='none',
                margins=c(8,20),
                symbreaks=FALSE,
                symm=FALSE,
                cex.main=0.75,
                density.info="none"
                )
      dev.off()
    },contentType='image/pdf')
  
  hide(id = "loading_page", anim = TRUE, animType = "fade")    
  show("mainTabsetPanel")
}

########################################################################################
#  Setaria Shiny Application - run app.
########################################################################################

# Run the application 
shinyApp(ui = ui, server = server)

########################################################################################
