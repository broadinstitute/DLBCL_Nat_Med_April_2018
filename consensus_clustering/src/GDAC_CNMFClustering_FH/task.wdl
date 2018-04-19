composite_task GDAC_CNMFClustering_FH {
  step GDAC_TopgenesforCluster[version=39] {
    output: File("${outputprefix}.expclu.gct") as expresion_file;
  }
  step GDAC_NmfConsensusClustering[version=49] {
    input: expfile=expresion_file;
    output: File("${outputprefix}.consensus.plot.k2.png") as plotk2, File("${outputprefix}.consensus.plot.k3.png") as plotk3, File("${outputprefix}.consensus.plot.k4.png") as plotk4, File("${outputprefix}.consensus.plot.k5.png") as plotk5, File("${outputprefix}.consensus.plot.k6.png") as plotk6, File("${outputprefix}.consensus.plot.k7.png") as plotk7, File("${outputprefix}.consensus.plot.k8.png") as plotk8, File("${outputprefix}.cophenetic.coefficient.txt") as coefficientfile, File("${outputprefix}.membership.txt") as membershipfile;
  }
  step GDAC_CNMFselectcluster[version=54] {
    input: inputexp=expresion_file, clumembership=membershipfile, cophenetic=coefficientfile;
    output: File("${output}.silfig.png") as sigkclus, File("${output}.bestclus.txt") as bestcluster, File("${output}.subclassmarkers.txt") as markerfile, File("${output}.cormatrix.png") as cormatixpng, File("${output}.selectmarker.txt") as selectmarkers, File("${output}.geneheatmap.png") as geneheatmap, File("${output}.geneheatmaptopgenes.png") as geneheatmapall;
  }
  step GDAC_CnmfReports[version=34] {
    input: expdata=expresion_file, kclus=sigkclus, bestclu=bestcluster, allcluster=membershipfile, markers=markerfile, cormatrix=cormatixpng, file.gif.2=plotk2, file.gif.3=plotk3, file.gif.4=plotk4, file.gif.5=plotk5, file.gif.6=plotk6, file.gif.7=plotk7, file.gif.8=plotk8, markersP=selectmarkers, heatmap=geneheatmap, heatmapall=geneheatmapall;
  }
}