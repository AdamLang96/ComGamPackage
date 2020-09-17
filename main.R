## still needs EB option
## could include nonlinear adjustment option

ComGamHarm <- function(feature.data, 
                       covar.data, 
                       training.indicies = NULL,
                       smooth.terms = NULL,
                       k.val = NULL,
                       verbose = TRUE,
                       model.diagnostics = FALSE) {
  
  ### check object types and parameters ###
  #include check for STUDY variable and include class check for factor
  if(!is.data.frame(feature.data)) {
    stop('ComGamHarm: feature.data is not of type "data.frame"')
  }
  if(!is.data.frame(covar.data)) {
    stop('ComGamHarm: covar.data is not of type "data.frame"')
  }
  if(!("STUDY" %in% colnames(covar.data))) {
    stop('ComGamHarm: covar.data must include column "STUDY"')
  }
  if(!is.factor(covar.data[["STUDY"]])) {
    stop('ComGamHarm: "STUDY" column must be of class "factor"')
  }
  if(any(is.na(feature.data)) | any(is.na(covar.data))) {
    stop('ComGamHarm: data frames cannot contain missing data')
  }
  if(!(nrow(feature.data) == nrow(covar.data))) {
    stop('ComGamHarm: # of rows inconsistent between data frames')
  }
  if(!is.null(smooth.terms) | !is.null(k.val)) {
    if(!(!is.null(smooth.terms) & !is.null(k.val))) {
      stop('ComGamHarm: both smooth.terms & k.val must be supplied if one is supplied')
     }
    if(!(length(smooth.terms == length(k.val)))) {
      stop('ComGamHarm: smooth.terms & k.val must be vectors of same length')
     
    }
  }
  data.dict        <-  BuildDict(covar.data =  covar.data)
  
  model.formula    <-  BuildFormula(covar.data   = covar.data,
                                    smooth.terms = smooth.terms,
                                    k.val        = k.val)
  
  models.list      <-  FitModel(feature.data     = feature.data,
                                covar.data       = covar.data,
                                model.formula    = model.formula,
                                training.indices = training.indicies,
                                verbose = verbose)
  
  stan.dict        <-  StanAcrossFeatures(feature.data = feature.data,
                                          covar.data   = covar.data,
                                          models.list  = models.list,
                                          data.dict    = data.dict)
  
  features.adj     <-  CalcGammaDelta(stan.dict =  stan.dict,
                                      data.dict = data.dict)
  
  features.results <-  ApplyGammaDelta(stan.dict   = stan.dict,
                                       site.params = features.adj,
                                       data.dict   = data.dict)
  
  if(model.diagnostics) {
    
    mod.diags <- ModelDiagnostics(mod.list = models.list)
    
    return.list <- list("harm.results"       = features.results,
                        "stan.dict"          = stan.dict,
                        "shift.scale.params" = features.adj,
                        "models.list"        = models.list,
                        "model.formula"      = model.formula,
                        "data.dict"          = data.dict,
                        "model.diagnostics"  = mod.diags)
  } else {
    
    return.list <- list("harm.results"       = features.results,
                        "stan.dict"          = stan.dict,
                        "shift.scale.params" = features.adj,
                        "models.list"        = models.list,
                        "model.formula"      = model.formula,
                        "data.dict"          = data.dict)
  }
  
  return(return.list)
}



