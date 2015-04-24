# prompted wrapper for Rotopol.qPLM.R, meant to be run as a batch file from Windows

Rotopol.tofile<-function(filename="test.csv"){
  source('~/GitHub/microtransit/Input-and-Imaging/Rotopol.qPLM.R')
  source('~/GitHub/microtransit/Analysis/qPLMtabulate.R')
  sample.name<-readline("Name the sample:")
  thickness<-as.numeric(readline("Specimen thickness (in microns):"))
  wavelength<-as.numeric(readline("Illumination wavelength (in nm--hit ENTER for default of 532):"))
  if (is.na(wavelength)) {
    wavelength<-532
  }
  birefringence<-as.numeric(readline("Birefringence of tissue (hit ENTER for extant bone default of 0.005):"))
  if (is.na(birefringence)) {
    birefringence<-0.005
  }
  cal<-Rotopol.qPLM(sample.name=sample.name,thickness=thickness,wavelength=wavelength,birefringence=birefringence,mask=TRUE, pics.out=FALSE, comp.bkgrnd=100)
  result<-qPLMtabulate(cal)
  write.csv(result, file = filename,row.names=FALSE,col.names=FALSE)
}
