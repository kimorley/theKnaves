plotCluster <- function(data,model,save=F,print=T){
	testData <- data.frame(t=seq(0,1,0.01))
	testData <- cbind(testData,r=predict(model$model,testData))
	clusterPlot <- ggplot() + 
			layer(
					data = data, mapping = aes(x = t, y = r, colour = call, alpha = conf),
					geom = "point"
			) + 
			layer(
					data = testData, mapping = aes(x = t, y = r),
					geom = "line"
			) +
			scale_colour_brewer(palette="Set1") + ylab("R") + xlab(expression(theta))
	if (print){
		print(clusterPlot)
	}
	if (save){
		ggsave(file="clusterPlot.pdf",height=5,width=5)
	}
}