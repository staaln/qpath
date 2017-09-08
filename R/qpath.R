#' qpath function
#'
#' qpath function
#' @export
#' @import pathifier
#' @examples

qpath<-function(data, allgenes, syms, pathwaynames, normals = NULL, ranks = NULL, attempts = 100, maximize_stability = TRUE, logfile = "", samplings = NULL, min_exp = 4, min_std = 0.4){
	PDS<-quantify_pathways_deregulation(data, allgenes, syms, pathwaynames, normals, ranks, attempts, maximize_stability, logfile, samplings, min_exp, min_std)
	p<-length(PDS$scores)
	qpathscores<-list(length=p)
	for (i in 1:p){
		print(i)
		a<-unlist(PDS$scores[[i]])
		m<-which.min(a)
		qpathscores[[i]]<-as.matrix(dist(PDS$curves[[i]][,PDS$compin[[i]]]))[m,]/max(as.matrix(dist(PDS$curves[[i]][,PDS$compin[[i]]]))[m,])
	}	
	names(qpathscores)<-pathwaynames
	list(qpathscores=qpathscores,PDS=PDS)
}