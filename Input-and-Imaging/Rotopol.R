# prompted wrapper for Rotopol.qPLM.R, meant to be run as a batch file from Windows

Rotopol<-function(){
  sample.name<-readline("Name the sample:")
  thickness<-as.numeric(readline("Specimen thickness (in microns):"))
  wavelength<-as.numeric(readline("Illumination wavelength (in nm--hit ENTER for default of 532):"))
  if (is.na(wavelength)) {
    wavelength<-532
  }
  birefringence<-as.numeric(readline("Birefringence of tissue (hit ENTER for Lab Turkey Tendon default of 0.0025):"))
  if (is.na(birefringence)) {
    birefringence<-0.0025
  }
  Rotopol.qPLM(sample.name=sample.name,thickness=thickness,wavelength=wavelength,birefringence=birefringence,mask=TRUE, comp.bkgrnd=100)
}
