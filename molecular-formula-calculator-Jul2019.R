##########################################
# NRC Molecular formula calculator 
# Author: Juris Meija, NRC Canada (2019)
# version: v.1.01 (30 July 2019)

options(repos = BiocInstaller::biocinstallRepos())
getOption("repos")

require(ecipex)
require(stringr)
require(Rdisop)
require(CHNOSZ)
require(shiny)
require(shinyjs)
require(rhandsontable)
require(DT)

#Isotope abundance dataframe (TICE-2013)
x.data <- read.csv("data/tice-2013-dataframe.csv")
data(nistiso)

# Helper functions to extract (1) nuclide masses, (2) isotopic abundances and (3) abundance covariances
mass <- function(el)   {x.data[which(x.data$element==el),]$mass}
x.mean <- function(el) {x.data[which(x.data$element==el),]$abundance}

# consolidate isotopologues that are within the specified mass window
compact = function(z, tol=0.1){
  mask = rep(0, length(z$mass))
  k = 0
  for (i in 1:length(z$mass)) { k = k+1; mask[abs(z$mass[i]-z$mass) <= 2*tol] = k }
  d = data.frame(cbind('mask'=mask, 'm'=z$mass, 'x'=z$abundance))
  # abundance-weighted mean for mass
  m.aggregate = by(d, d$mask, function(x) weighted.mean(x$m, x$x))
  # sum isotopic abundances
  x.aggregate = by(d, d$mask, function(x) sum(x$x))
  list(m=m.aggregate, x=x.aggregate)
}

mol.formula = function(x, constraint='CHNO', elements, enr = TRUE, sen = TRUE, tol.ppm = 3, df.enr, charge){
  # parse input data
  m.nat=x$m.nat
  x.nat=x$x.nat
  m.enr=x$m.enr
  x.enr=x$x.enr
  charge=as.double(charge)
  # helper functions ##################################################
  # consolidate all isotopologues that are within the specified mass window
  compact = function(z, tol=0.1){
    mask = rep(0, length(z$mass))
    k = 0
    for (i in 1:length(z$mass)) { k = k+1; mask[abs(z$mass[i]-z$mass) <= 2*tol] = k }
    d = data.frame(cbind('mask'=mask, 'm'=z$mass, 'x'=z$abundance))
    # abundance-weighted mean for mass
    m.aggregate = by(d, d$mask, function(x) weighted.mean(x$m, x$x))
    # sum isotopic abundances
    x.aggregate = by(d, d$mask, function(x) sum(x$x))
    list(m=m.aggregate, x=x.aggregate)
  }
  # check if the isotopic pattern of the natural substance and its hypothesis are in agreement (R^2 >0.9) 
  goodness.nat = function(formula){
    #if(is.na(formula))return(0)
    z=ecipex(formula, isoinfo = x.data, sortby = "mass", limit=1e-4)[[1]]
    d=compact(z)
    m.j = outer(d$m, m.nat, function(x, y) abs(x-y))
    matches = which(m.j < 0.1, arr.ind=TRUE)
    # if(length(matches)==0|any(is.na(matches))) return(0)
    result=cor(d$x[matches[,1]], x.nat[matches[,2]])
    names(result)=formula
    return(result)
  }
  # check if the isotopic pattern of the enriched substance and its hypothesis are in agreement (R^2 >0.9) 
  goodness.enr = function(formula){
    #if(is.na(formula)) return(result=0)
    z=ecipex(formula, isoinfo = df.enr, sortby = "mass", limit = 1e-4)[[1]]
    d=compact(z)
    m.j = outer(d$m, m.enr, function(x, y) abs(x-y))
    matches = which(m.j < 0.1, arr.ind=TRUE)
    # if(length(matches)==0|any(is.na(matches)))return(FALSE)
    result=cor(d$x[matches[,1]], x.enr[matches[,2]])
    names(result)=formula
    return(result)
  }
  # match molecular formulas against the specified constraints
  constraint.f = function(formula, constraint) {
    constraint.split = str_extract_all(constraint,"([A-Z][a-z]*)([0-9]*)")[[1]]
    p1 = '(?:[A-Z]|$)'
    p2 = '(?:[A-Z]|[0-9]|$)'
    all(sapply(constraint.split, function(p) grepl(pattern=paste0(p,ifelse(grepl('[0-9]',p),p1,p2)), x=formula)))
  }
  
  # convert multi-charges to zero charge
  m.nat = m.nat*abs(charge) + 0.0005486*charge
  m.enr = m.enr*abs(charge) + 0.0005486*charge
  
  # step 1
  e = str_extract_all(elements, "([A-Z][a-z]*)")[[1]]
  qq=q=decomposeIsotopes(m.nat, x.nat, 
                         elements = initializeElements(e), 
                         ppm=tol.ppm, mzabs = 0, z=0, maxisotopes = 10)
  n.1 = length(q$formula)
  
  # step 2
  # eliminate formulas that do not match the provided constraints
  #if(constraint=='') formulas=q$formula[grepl("N",q$formula)&grepl("C",q$formula)&grepl("H",q$formula)&grepl("O",q$formula)]
  #if(constraint!='') formulas=q$formula[grepl(constraint,q$formula)&grepl("N",q$formula)&grepl("C",q$formula)&grepl("H",q$formula)&grepl("O",q$formula)]
  formulas=q$formula[sapply(q$formula, constraint.f, constraint=constraint)]
  n.2 = length(formulas)
  
  if(n.2<1)return(list(formula=NULL, score=NULL, cor.N = NULL, cor.E = NULL, n=list(NULL,NULL,NULL,NULL,NULL,NULL)))  
  # step 3
  # implement Senior rules for neutral formulas
  # 1 function to add or subtract the hydrogen
  f.modify = function(x,charge) {
    f = makeup(x)
    f[names(f)=="H"] = f[names(f)=="H"] - charge
    paste0(names(f), f, collapse="")
  }
  
  # 2 Senior rule check, doi: 10.2307/2372318 and http://khimiya.org/volume12/vol12_6_morikawa.pdf
  ### 1 sum of valences is even OR the total number of atomc having odd valence is even
  ### 2 sum of valences is >= twice the number of atoms minus 1
  senior = function(f){
    x=makeup(f)
    nC = ifelse(length(x[names(makeup(x))=="C"])==0,0,x[names(makeup(x))=="C"])
    nH = ifelse(length(x[names(makeup(x))=="H"])==0,0,x[names(makeup(x))=="H"])
    nN = ifelse(length(x[names(makeup(x))=="N"])==0,0,x[names(makeup(x))=="N"])
    nO = ifelse(length(x[names(makeup(x))=="O"])==0,0,x[names(makeup(x))=="O"])
    nS = ifelse(length(x[names(makeup(x))=="S"])==0,0,x[names(makeup(x))=="S"])
    rule1 = unname(((4*nC+1*nH+3*nN+2*nO+2*nS) %%2 == 0) | ((nH + nN) %%2 == 0))
    rule2 = unname((4*nC+1*nH+3*nN+2*nO+2*nS) >= (2*(nC + nH + nN + nO + nS) - 1))
    rule1 & rule2
  }
  n.3 = n.2
  if(sen){
    formulas.neutral = unname(sapply(formulas, f.modify, charge=charge))
    formulas.final = formulas[sapply(formulas.neutral, senior)]
    formulas <- formulas.final
    n.3 = length(formulas)
  }
  
  if(n.3<1)return(list(formula=NULL, score=NULL, cor.N = NULL, cor.E = NULL, n=list(NULL,NULL,NULL,NULL,NULL,NULL)))
  # step 4
  # eliminate formulas whose natural isotope patterns disagree
  cor.nat = sapply(formulas, goodness.nat)
  formulas <- formulas[cor.nat > 0.9]
  n.4 = length(formulas)
  
  # steps 5-6
  n.5=n.6=n.4
  
  if(enr==TRUE){
    # eliminate formulas whose labeled isotope patterns disagree
    cor.enr = sapply(formulas, goodness.enr)
    formulas <- formulas[cor.enr > 0.9]
    n.5 = length(formulas)

    # eliminate formulas whose masses are further than +/- k*tol.ppm away for the enriched substance
    m.enr.correct = sapply(formulas, function(f, k=2){
      if(is.na(f)) return(FALSE)
      m.enr.theor = ecipex(f, isoinfo=df.enr, sortby="abundance", limit = 1e-4)[[1]]$mass
      mztol = mean(abs(m.enr.theor*k*tol.ppm*1e-6))
      all(sapply(m.enr, function(x) any(abs(x - m.enr.theor) < mztol)))
    })
    
    formulas <- formulas[m.enr.correct]
    n.6 = length(formulas)
  }
  
  f = na.omit(formulas)
  if(length(f)<1)return(list(formula=NULL, score=NULL, cor.N = NULL, cor.E = NULL, n=list(NULL,NULL,NULL,NULL,NULL,NULL)))
  
  # average deviation in masses
  d.nat=d.enr=c()
  for(qqq in f){
    a0=ecipex(qqq,isoinfo=x.data,sortby="mass",limit=1e-5)[[1]]
    a1=ecipex(qqq,isoinfo=df.enr,sortby="mass",limit=1e-5)[[1]]
    # theoretical spectra
    pat.nat = compact(a0)
    pat.enr = compact(a1)
    dm.nat=dm.enr=c()
    for(i in 1:length(m.nat)) dm.nat[i]=min(abs(m.nat[i] - pat.nat$m))
    for(i in 1:length(m.enr)) dm.enr[i]=min(abs(m.enr[i] - pat.enr$m))
    dm.1 = 1000*mean(dm.nat)
    dm.2 = 1000*mean(dm.enr)
    d.nat = c(d.nat, dm.1)
    d.enr = c(d.enr, dm.2)
  }

  ### Report the results
  d.nat = formatC(d.nat,format='f',digits=2)
  d.enr = ifelse(enr,formatC(d.enr,format='f',digits=2), NA )
  z = qq$score[qq$formula %in% f]
  # calculate the relative score
  score=formatC(z/max(z, na.rm=TRUE),format='f',digits=3)
  #score=formatC(z,format='f',digits=4)
  cor.N = formatC(sapply(f, goodness.nat),format='f',digits=3)
  cor.E = ifelse(enr,formatC(sapply(f, goodness.enr),format='f',digits=3), NA)
  return(list(formula=f, score=score, cor.N = cor.N, cor.E = cor.E, dm.N = d.nat, dm.enr = d.enr, n=list(n.1,n.2,n.3,n.4,n.5,n.6)))
}

# plot the results
f.plot = function(f,score,m.nat,x.nat,m.enr,x.enr,min=NULL,max=NULL,df.N,df.E,enr,charge){
  if(length(f)<1)return(NULL)
  require(gsubfn)
  a0=ecipex(f,isoinfo=df.N,sortby="mass",limit=1e-5)[[1]]
  a1=ecipex(f,isoinfo=df.E,sortby="mass",limit=1e-5)[[1]]
  # theoretical spectra
  pat.nat = compact(a0)
  pat.enr = compact(a1)
  
  pat.nat$m <- (pat.nat$m - 0.0005486*charge)/abs(charge)
  pat.enr$m <- (pat.enr$m - 0.0005486*charge)/abs(charge)
    
  plot(pat.nat$m,pat.nat$x,t="h",yaxt="n",xlim=c(ifelse(is.null(min),min(pat.nat$m)-1,min),ifelse(is.null(max),max(pat.enr$m)+1,max)), main="", xlab="Atomic mass, Da", ylab="", cex.main=0.9)
  #mtext(parse(text=substr(ss,1,nchar(ss)-1)), side=3, line=0.3, cex=1, font=2, adj=0)
  mtext(f, side=3, line=0.3, cex=1, font=2, adj=0)
  mtext(paste('score:',score), side=3, line=0.3, cex=0.9, font=2, adj=1)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray95")
  if(enr)segments(pat.enr$m,0,pat.enr$m,pat.enr$x, col="#C93312",lwd=2.5)
  segments(pat.nat$m,0,pat.nat$m,pat.nat$x,col=1,lwd=2.5)
  # experimental spectra
  if(enr)points(m.enr, x.enr*max(pat.enr$x)/max(x.enr), pch=19, col='#C93312', cex=2) # labeled
  points(m.nat, x.nat*max(pat.nat$x)/max(x.nat), pch=19, col=1, cex=2) # natural
  box()
}

mass.deviation = function(f,m.nat,x.nat,m.enr,x.enr,df.N,df.E,enr){
  # average deviation in masses
  a0=ecipex(f,isoinfo=df.N,sortby="mass",limit=1e-5)[[1]]
  a1=ecipex(f,isoinfo=df.E,sortby="mass",limit=1e-5)[[1]]
  # theoretical spectra
  pat.nat = compact(a0)
  pat.enr = compact(a1)
  
  dm.nat=dm.enr=c()
  for(i in 1:length(m.nat)) dm.nat[i]=min(abs(m.nat[i] - pat.nat$m))
  for(i in 1:length(m.enr)) dm.enr[i]=min(abs(m.enr[i] - pat.enr$m))
  dm.1 = 1000*mean(dm.nat)
  dm.2 = 1000*mean(dm.enr)
  if(enr)return(as.data.frame(list(dm.1,dm.2)))
  return(as.data.frame(list(dm.1,0)))
}


## TEST DATA of microcystin [Asp3]MC-RCit (C48H71N12O13) with x(15N) = 0.99 (z=-1)
   init.df.N = data.frame(mass=c(1023.5287, 1024.5318, 1025.5349, 1026.5378),
                       abundance=c(558081, 305344, 86862, 17989))
   init.df.E = data.frame(mass=c(1034.4960, 1035.4934, 1036.4963, 1037.4992),
                          abundance=c(110920, 860953, 447304, 130954))

### server file
server <- function(input, output, session) {
  
  output$A <- renderText(paste0("Isotopic abundance of ", c('carbon', 'nitrogen')[c(input$element=='C', input$element=='N')],"-",c(13, 15)[c(input$element=='C', input$element=='N')]))

  dfiso <- x.data

  # hot data input table (mol formula calculator natural)
  output$hot.N <- renderRHandsontable({
    if (is.null(input$hot.N)) { DF = init.df.N } else { DF = hot_to_r(input$hot.N) }
    rhandsontable(DF, readOnly = FALSE, stretchH = "all", selectCallback = TRUE) %>%
      hot_context_menu(allowColEdit = FALSE ) %>%
      hot_validate_numeric(cols = 1:2, min = 0) %>%
      hot_cols('float', format = '0.000')
  })
  # hot data input table (mol formula calculator labeled)
  output$hot.E <- renderRHandsontable({
    if (is.null(input$hot.E)) { DF = init.df.E } else { DF = hot_to_r(input$hot.E) }
    rhandsontable(DF, readOnly = FALSE, stretchH = "all", selectCallback = TRUE) %>%
      hot_context_menu(allowColEdit = FALSE ) %>%
      hot_validate_numeric(cols = 1:2, min = 0) %>%
      hot_cols('float', format = '0.000')
  })
  
  user.data  <- reactive({
    if (is.null(input$hot.N)&is.null(input$hot.E)) { z.N = init.df.N; z.E = init.df.E } else { z.N = hot_to_r(input$hot.N); z.E = hot_to_r(input$hot.E)}
    list(m.nat=z.N[,1], x.nat=z.N[,2], m.enr=z.E[,1], x.enr=z.E[,2], tol.ppm = input$tolerance)
  })
  
  
  observeEvent(input$button, {
            
            df0 = x.data
            df0$abundance[df0$element==input$element] <- c(1-input$abundance, input$abundance)
            
            res=mol.formula(user.data(), enr=ifelse(input$enr=='Yes',TRUE,FALSE), sen=ifelse(input$sen=='Yes',TRUE,FALSE), elements=input$elements,constraint = input$constraints, tol.ppm = input$tolerance, df.enr = df0, charge=as.double(paste(input$charge)))
            res1 = as.data.frame(res)
            # error : no matching peaks
            if(length(res1$formula)<1){
              output$text3  <- renderText({ paste("<font color=\"#FF0000\">", "NO MATCHING SUBSTANCES", "</font>") })
              output$text4  <- renderText({ paste("<font color=\"#FF0000\">", "Verify your input data or try increasing the m/z tolerance", "</font>") })
            }
            
            if(length(res1$formula)>0){
            
            output$text1  <- renderText({paste("Molecular formulas that match the experimental data and all computational constraints") })
            output$T1  <- DT::renderDataTable(res1[,c(1:6)], selection = list(mode='single', selected=1), rownames= FALSE, server = FALSE)
            output$text3  <- renderText({ paste("Molecular formula matches: ","(step 1)",res1[1,7],"(step 2)",res1[1,8],"(step 3)",res1[1,9],"(step 4)",res1[1,10],"(step 5)",res1[1,11],"(step 6)",res1[1,12]," (see notes for explanation)") })
            output$text5  <- renderText({ paste("Click on the molecular formula to see the experimental/calculated isotopic patterns.") })
          
            output$plot <- renderPlot({
              s = input$T1_rows_selected
              if (length(s)) {
                f.plot(paste(res1[s,1]),paste(res1[s,2]),
                       user.data()$m.nat,user.data()$x.nat,user.data()$m.enr,user.data()$x.enr,
                       df.N=x.data,df.E=df0,enr=ifelse(input$enr=='Yes',TRUE,FALSE),charge=as.double(paste(input$charge)))
                }
            })
            
            
            #output$T2 <- DT::renderDataTable({
            #    s = input$T1_rows_selected
            #    if(length(s))
            #    mass.deviation(paste(res1[s,1]),user.data()$m.nat,user.data()$x.nat,user.data()$m.enr,user.data()$x.enr,
            #                   df.N=x.data,df.E=df0,enr=ifelse(input$enr=='Yes',TRUE,FALSE))
            #  }, selection = 'single', rownames= FALSE, server = FALSE)
            
            
            }
            
            })
  

  shinyjs::onclick("toggleextra.mol", shinyjs::toggle(id = "filterextra.mol", anim = TRUE))
  shinyjs::onclick("togglenotes.mol", shinyjs::toggle(id = "filternotes.mol", anim = TRUE))
  
}

### user interface
ui <- fluidPage(
  #tags$head(tags$link(rel="shortcut icon", href="URL-to-favicon"))
  titlePanel( title="Molecular formula calculator" ),
  shinyjs::useShinyjs(),
  sidebarLayout(
    
    sidebarPanel(
      #column(6,selectInput("charge", label = "Select the charge-number", choices=seq(-3,+3), selected=-1, width = '100%')),
      #column(6,textInput("adduct", label = "Enter the adduct", value="H", width = '100%')),
      selectInput("charge", label = "Select the charge-state", choices=c('-2 [M-2H]2-'=-2,'-1 [M-H]-'=-1,'+1 [M+H]+'=1,'+2 [M+2H]2+'=2), selected=-1, width = '85%'),
      h5(tags$b("The observed mass spectrum of the natural substance")),
      rHandsontableOutput("hot.N"),
      helpText("right-click to add or delete rows"),
      br(),
      h5(tags$b("The observed mass spectrum of the labeled substance")),
      rHandsontableOutput("hot.E"),
      helpText("right-click to add or delete rows"),
      br(),
      selectInput("element", label = "Select the isotopically labeled element", choices = c('carbon'='C','nitrogen'='N'), selected = 'N', width = '85%'),
      sliderInput("abundance", label = textOutput("A"), min=0, max=1, value=0.99, round=TRUE, step=0.01, width='85%'),
      h5("Additional parameters ", a(id = "toggleextra.mol", "show/hide")),
      shinyjs::hidden(div(id = "filterextra.mol",
                          fluidRow(
                            column(12, numericInput("tolerance", label = "Enter m/z tolerance (ppm)", value = 3, max = 50, min=1, step = 1, width = '85%')),
                            column(12, textInput("elements", "What elements to consider", value="CHNOS", placeholder = "for example: CHNO or CHNOS", width = '85%')),
                            column(12, textInput("constraints", "Constraints to the molecular formula", value="CHNO", placeholder = "for example: CHNO or N12C48 or C48N12", width = '85%')),
                            column(12, radioButtons("enr", "Should the isotope pattern of labeled substance be considered at all?", choices=c("Yes","No"), selected="Yes", inline = TRUE, width = '85%')),
                            column(12, radioButtons("sen", "Should Senior rules be applied?", choices=c("Yes","No"), selected="Yes", inline = TRUE, width = '85%'))
                          ))),
      br(),
      br(),
      conditionalPanel(condition = "!$('html').hasClass('shiny-busy')",   
                         actionButton("button", label = "Perform data analysis", icon = icon('bar-chart-o'))),
      conditionalPanel(condition = "$('html').hasClass('shiny-busy')",   
                         actionButton("button", label = "busy...", icon = icon('hourglass-half'))),
      br(),
      p("NRC Molecular Formula Calculator (2019) v.1.01")
    ),
    
    mainPanel(
      tags$style(type="text/css",
                 ".shiny-output-error { visibility: hidden; }",
                 ".shiny-output-error:before { visibility: hidden; }"),
      fluidRow(column(helpText("This calculator evaluates two experimental mass spectra of a given substance: one with natural and the other with isotopically altered isotopic composition. This information is used to determine the likely molecular formula."),width=11)),
      br(),
      br(),
      #conditionalPanel(condition = "!$('html').hasClass('shiny-busy')",
                       htmlOutput("text3"),br(),br(),
                       strong(htmlOutput("text1")),br(),br(),
                       fluidRow(column(DT::dataTableOutput("T1"),width=11)),
                       htmlOutput("text5"),
                       htmlOutput("text4"),
                       
      #           ),
      fluidRow(column(plotOutput(("plot")),width=11)),
      #fluidRow(column(DT::dataTableOutput("T2"),width=11)),
      br(),
      conditionalPanel(condition = "!$('html').hasClass('shiny-busy')",
      h5("Notes and explanations", a(id = "togglenotes.mol", "show/hide")),
      shinyjs::hidden(div(id = "filternotes.mol",
                          fluidRow(column(
                            p("This calculator is based mainly on R packages ecipex (fast calculation of isotope patterns using Fourier transform) and Rdisop (mass decomposition of isotope patterns using efficient money changing algorithm). Several other R packages were used in this calculator."),
                            p("Candidate molecular formulas are produced in six steps: (1) All matches for the natural MS data are produced using Rdisop package within the given mass tolerance constraints, (2) the candidates are reduced to eliminate formulas with no C, N, H, or O atoms, or to comply with the provided constrains, (3) candiate formulas are reduced to comply with the Senior rules of molecular composition, (4) formulas whose natural isotope patterns disagree with the expected model (ecipex) are discarded, (5) formulas whose labeled isotope patterns disagree with the expected model (ecipex) are discarded, and (6) candidates whose masses in the labeled spectra deviate from the expected model as dictated by the known isotopic abundance of the labeled medium are discarded."),
                            p("Created using R and Shiny by Juris Meija (2019) NRC Canada"),width=11)
                          )))
      )
      
    )
  )
)

shinyApp(ui = ui, server = server)
