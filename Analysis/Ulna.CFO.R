# 21 FEB '15 bird CFO analysis function
# pull all rotopol into qPLM objects first

Ulna.CFO<-function(sample.name) {
  .data<-qPLMtabulate(sample.name)
  .mech.mod<-gen.mech.model(.data)
  .MEM<-gen.sp.MEM(.data,cortical=TRUE,maps=min(c(nrow(.data$distance)-2,300)))
  .sp.mod<-dbMEM.model.scalogram(.MEM)
  parts<-list()
  parts<-c(".data",".mech.mod",".MEM",".sp.mod")
  save(list=parts,file=paste(sample.name,".ulna.R",sep=""))
}