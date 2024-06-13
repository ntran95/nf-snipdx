classifySample <- function(sample.vec, cluster.centers, cluser.df, sample_type, plot=NULL, verbose=TRUE) {

    # sample_type: blood or FFPE



    v.norm <- sample.vec/sum(sample.vec, na.rm=TRUE)
    v.norm[is.na(v.norm)] <- 0
    euclidean <- function(a, b) sqrt(sum((a - b)^2))
    cluster.dists <- c()
    for (i in 1:nrow(cluster.centers)) {
        cluster.dists[i] <-  euclidean(v.norm, cluster.centers[i,])
    }



    names(cluster.dists) <- rownames(cluster.centers)




    if (sample_type=='FFPE' | sample_type=='ffpe') {
        cluster.dists <- cluster.dists[grepl('ffpe',names(cluster.dists))]
    }  else if (sample_type=='Blood' | sample_type=='blood' | sample_type=='BLOOD') {
        cluster.dists <- cluster.dists[grepl('blood',names(cluster.dists))]
    }

    cluster.assignment.no <- which.min(cluster.dists)
    cluster.assignment.ascCode <- names(cluster.dists)[cluster.assignment.no]
    if (verbose) {
        print('cluster distances')
        print(sort(cluster.dists))
        print(paste('classification', cluster.assignment.ascCode))
   }

    return(cluster.assignment.ascCode)

}