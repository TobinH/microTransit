# 21 FEB '15 bird CFO analysis function
# pull all rotopol into qPLM objects first

Ulna.CFO<-function(sample.name) {
  .data<-qPLMtabulate(sample.name)
  .scatter<-scatter.matrix.summary(.data)
  .mech.mod<-gen.mech.model(.data)
  .MEM<-gen.sp.MEM(.data,cortical=TRUE,maps=min(c(nrow(.data$distance)-2,300)))
  .sp.mod<-dbMEM.model.scalogram(.MEM)
  results<-list(data=.data,scatter=.scatter,mech.mod=.mech.mod,MEM=.MEM,sp.mod=.sp.mod)
  return(results)
}