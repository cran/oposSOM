
run.pipeline = function( indata=NULL, dataset.name = "Unnamed", dim.som1 = 20, dim.som2 = 20 )
{

	preferences = list()
	preferences$dataset.name = dataset.name
	preferences$dim.som1 = dim.som1
	preferences$dim.som2 = dim.som2
			
	preferences$sample.spot.cutoff = 0.65
	preferences$summary.spot.cutoff = 0.95
		
	



	cat( "\nStarted:", format(Sys.time(), "%a %b %d %X\n" ) )
	cat( "Setting:", preferences$dataset.name, "\n" ) 
	cat( "1SOM Dim:", preferences$dim.som1,"\n2SOM Dim:",preferences$dim.som2,"\n\n" )
	flush.console()




	## Preprocessing input data object


	if( class( indata ) != "matrix" || mode( indata ) != "numeric" )
	{
		stop("Input data must be numerical matrix.")
	}

	if( length( rownames( indata ) ) == 0 )
	{
		rownames( indata ) = as.character( 1:nrow(indata) )
	}
	
	if( length( colnames( indata ) ) == 0 )
	{
		colnames( indata ) = paste( "Sample", c( 1:ncol(indata) ) )
	}
                  
	na.rows = which( apply( apply( indata, 1, is.na ), 2, sum ) > 0 )
	if( length( na.rows ) > 0 )
	{
		indata = indata[ -na.rows, ]
		cat( "\n!!Removed NAs from data set!!\n" ); flush.console()
	}	

	indata = indata - apply( indata, 1, mean )






	## SOM Processing

	cat( "\nProcessing SOM\n" );	flush.console()

	som.result = som( indata, xdim=preferences$dim.som1, ydim=preferences$dim.som1 )


	metadata = som.result$code
	colnames(metadata) = colnames( indata )

	som.result$code = NA






	## set up global variables

	group.labels = rep( "sample", ncol( indata ) )
	names(group.labels) = colnames( indata )

	group.colors = rainbow( ncol( indata )  )
	names(group.colors) = colnames( indata )

	colnames( indata ) = make.unique( colnames( indata ) )
	names( group.labels ) = make.unique( names( group.labels ) )
	names( group.colors ) = make.unique( names( group.colors ) )


	LETTERS = c( LETTERS,   as.vector( sapply( 1:10, function(x){ paste( LETTERS, x, sep=""  ) } ) )   )


	colramp = colorRampPalette( c( "darkblue","blue","lightblue","green","yellow","red","darkred" ) )


	som.nodes = ( som.result$visual[,"x"] + 1 ) +  som.result$visual[,"y"] * preferences$dim.som1
	names( som.nodes ) = rownames( indata )



	## load gene set annotations from ensembl database

	cat( "Load Annotation Data\n\n" ); flush.console()

	
#	if( length( .find.package( "biomaRt", quiet=T ) ) == 0 )
#	{
#		source("http://bioconductor.org/biocLite.R")
#		biocLite( "biomaRt" )
#	}
#	library( biomaRt )
	
	
	gene.names = rownames(indata)
	names(gene.names) = rownames(indata)
	gene.descriptions = rep("", nrow(indata))
	names(gene.descriptions) = rownames(indata)
	gene.ids = rownames(indata)
	names(gene.ids) = rownames(indata)
	gs.def.list = list()
	preferences$geneset.analysis = F
	
	
	mart = useMart('ensembl')
	ds = listDatasets(mart)$dataset
	use.dataset = ""
	if( length( grep("ENS", rownames(indata)) ) > 0 )
		for( i in ds )
		{
			try({
				mart<-useDataset( as.character(i),mart=mart)
	
				biomart.table = getBM( "ensembl_gene_id", "ensembl_gene_id", rownames(indata)[seq(1,nrow(indata),length.out=10)], mart, checkFilters=F )		
				if( nrow(biomart.table) != 0 )
				{
					use.dataset = as.character(i)		
					break
				}
			}, silent=T )
		}
	if( use.dataset == "" )
	{		
		cat("Could not download GO annotation. Gene set analysis will not be performed.\n\n")		
	} else
	{
		unique.protein.ids = unique(gene.ids)
		
		biomart.table = getBM( c( "ensembl_gene_id", "go_id", "name_1006", "namespace_1003" ) , "ensembl_gene_id", unique.protein.ids, mart, checkFilters=F )
		gs.def.list = tapply( biomart.table[,1], biomart.table[,3], c  )
		gs.def.list = lapply( gs.def.list, function(x){ list( Genes=x, Type="" ) } )
		
		gs.def.list = gs.def.list[ - which( names(gs.def.list) == "" ) ]
		gs.def.list = gs.def.list[ - which( sapply( sapply( gs.def.list, head, 1 ), length ) < 10 ) ]	
		
		gs.def.list = gs.def.list[ order( names( gs.def.list ) ) ]
		
		gs.def.list = lapply( gs.def.list, function(x)
		{ 
			x$Genes = intersect( x$Genes, unique.protein.ids )	
			return(x)
		} )
	
		empty.gs = which( sapply( sapply( gs.def.list, head, 1 ), length ) == 0 )
		if( length(empty.gs) > 0 )
		{
			gs.def.list = gs.def.list[ - empty.gs ]
		}
	
		if( length(gs.def.list) > 0 )
		{                                             	
			preferences$geneset.analysis = T
			cat( "Download of", length(gs.def.list), "GO sets using",use.dataset,"\n\n" ); flush.console()
			
		} else
		{
			cat("Could not download GO annotation. Gene set analysis will not be performed.\n\n")
		}		
		
	}
	




	## Output


	files.name = preferences$dataset.name


	while( file.exists( paste( files.name, ".html", sep="" ) ) )
		files.name = paste( files.name, "+", sep="" )


	dir.create( paste( files.name, "- Results" ), showWarnings=F )


	

		col.pix = function( m, x, y, col, threshold=NA )
		{
			hash.list = c( x+y*preferences$dim.som1 )
			pos.list = list( c(x,y) )
	
			while( length( hash.list ) > 0 )
			{
				x = pos.list[[1]][1]	
				y = pos.list[[1]][2]
	
				if( ( is.na(threshold) && !is.na(m[x,y]) && m[x,y] != col  ) || ( !is.na(threshold) && m[x,y] > threshold ) ) 
		      	{
		            	m[x,y] = col
	
			 		if( x-1 > 0 && length( intersect( (x-1)+(y)*preferences$dim.som1, hash.list ) ) == 0 )
					{
						hash.list = c( hash.list, (x-1)+(y)*preferences$dim.som1 )
						pos.list[[ length(pos.list)+1 ]] = c( x-1, y )
					}
			 		if( y-1 > 0 && length( intersect( (x)+(y-1)*preferences$dim.som1, hash.list ) ) == 0 )
					{
						hash.list = c( hash.list, (x)+(y-1)*preferences$dim.som1 )
						pos.list[[ length(pos.list)+1 ]] = c( x, y-1 )
					}
			 		if( x+1 <= nrow(m) && length( intersect( (x+1)+(y)*preferences$dim.som1, hash.list ) ) == 0 )
					{
						hash.list = c( hash.list, (x+1)+(y)*preferences$dim.som1 )
						pos.list[[ length(pos.list)+1 ]] = c( x+1, y )
					}
			 		if( y+1 <= ncol(m) && length( intersect( (x)+(y+1)*preferences$dim.som1, hash.list ) ) == 0 )
					{
						hash.list = c( hash.list, (x)+(y+1)*preferences$dim.som1 )
						pos.list[[ length(pos.list)+1 ]] = c( x, y+1 )
					}
	
				}
	
				hash.list = hash.list[ -1 ]
				pos.list = pos.list[ -1 ]
			}
	
			return( m )
	
		}	
	
	
	
	
		get.neighbors = function( x, y )
		{
			ret = list()
		
			if( x > 1 ) ret[['l']] = c( x-1, y )
			if( x < preferences$dim.som1 ) ret[['r']] = c( x+1, y )
			if( y > 1 ) ret[['u']] = c( x, y-1 )
			if( y < preferences$dim.som1 ) ret[['d']] = c( x, y+1 )
			
			return( ret )
		}



	GeneSet.HG = function( list.symbols, all.gene.symbols, gs.def.list, sort=F ) 
	{
	
		list.symbols = list.symbols[ !is.na( list.symbols ) ]
	
	
	
		N = length( all.gene.symbols )
		n = length( list.symbols )
	
	
		definition.p = rep( 0, length(gs.def.list) )
		names( definition.p ) = names(gs.def.list)
	
		for( i in 1:length( gs.def.list ) )
		{
			k = length( which( list.symbols %in% gs.def.list[[i]]$Genes ) )
			M = min( N, length( gs.def.list[[i]]$Genes ) )
	
			p = 0
			if( k > 0 )
				for( P.ki in 0:k )
					p = p + dhyper( P.ki, M, N-M, n )
			definition.p[ i ] = 1 - p
	
			if( definition.p[ i ] < 0 )		
				definition.p[ i ] = 0
		}	
	
	
		o = c(  1:length( definition.p )  )
		if( sort ) o = order( definition.p )
	
	
		ret = list()
	
		ret$p.value = definition.p[o]
	
		return( ret )
	
	}
	
	
	
	
		n.sets = 20
		n.samples = min( 20, ncol(indata) )
	
	
	
	
	
		##### Overexpression Cluster ######
	
	       
	
		
		peaks = apply( apply( metadata, 2, function(x){ ( x - min(x) ) / ( max(x) - min(x) ) } ), 1, max, na.rm=T )
	
	
	
		GS.infos.overexpression = list()
	
		GS.infos.overexpression$overview.map = peaks
	
	
	
		
		
	
	
		peaks = matrix( peaks, preferences$dim.som1, preferences$dim.som1 )
		peaks[ which( peaks < quantile( peaks, preferences$summary.spot.cutoff ) ) ] = NA
	
		e.cluster = matrix( NA, preferences$dim.som1, preferences$dim.som1 )
		e.cluster[ which( !is.na(peaks) ) ] = -1
		count.cluster = 1	
	
	
		while( any( !is.na( peaks ) ) )
		{
			
			start.pix = which( peaks == max(peaks,na.rm=T), arr.ind=T )[1,]
	
			e.cluster = col.pix( e.cluster, start.pix[1], start.pix[2], count.cluster )
	
			peaks[ which( e.cluster == count.cluster ) ] = NA
	
			count.cluster = count.cluster + 1
		}
	
	
	
		e.cluster = as.vector( e.cluster )
	
		for( i in 1:length( na.omit( unique( e.cluster ) ) ) )
		{
			nodes = which( e.cluster == i )
	
	
			geneset.genes = rownames( indata )[ which( som.nodes %in% nodes ) ]
			geneset.ids = unique( na.omit( gene.ids[ geneset.genes ] ) )
	
	
			GS.infos.overexpression[[ LETTERS[i] ]] = list()
	
			GS.infos.overexpression[[ LETTERS[i] ]]$metagenes = as.numeric( nodes )
			GS.infos.overexpression[[ LETTERS[i] ]]$genes = geneset.genes
	
			GS.infos.overexpression[[ LETTERS[i] ]]$mask = rep( NA, (preferences$dim.som1*preferences$dim.som1) )
			GS.infos.overexpression[[ LETTERS[i] ]]$mask[ as.numeric( nodes ) ] = 1
	
			if( preferences$geneset.analysis )
				GS.infos.overexpression[[ LETTERS[i] ]]$HG.p = GeneSet.HG( geneset.ids, unique.protein.ids, gs.def.list, sort=T )$p.value
		}
	
	
	
	
	
	
	
	
		##### Underexpression Cluster ######
	
	
	
		
		peaks = apply( apply( metadata, 2, function(x){ ( x - min(x) ) / ( max(x) - min(x) ) } ), 1, min, na.rm=T )
	
	
	
		GS.infos.underexpression = list()
	
		GS.infos.underexpression$overview.map = peaks
	
	
	
		
		
	
	
		peaks = matrix( peaks, preferences$dim.som1, preferences$dim.som1 )
		peaks[ which( peaks > quantile( peaks, 1-preferences$summary.spot.cutoff ) ) ] = NA
	
		e.cluster = matrix( NA, preferences$dim.som1, preferences$dim.som1 )
		e.cluster[ which( !is.na(peaks) ) ] = -1
		count.cluster = 1	
	
	
		while( any( !is.na( peaks ) ) )
		{
			
			start.pix = which( peaks == min(peaks,na.rm=T), arr.ind=T )[1,]
	
			e.cluster = col.pix( e.cluster, start.pix[1], start.pix[2], count.cluster )
	
			peaks[ which( e.cluster == count.cluster ) ] = NA
	
			count.cluster = count.cluster + 1
		}
	
	
	
		e.cluster = as.vector( e.cluster )
	
		for( i in 1:length( na.omit( unique( e.cluster ) ) ) )
		{
			nodes = which( e.cluster == i )
	
	
			geneset.genes = rownames( indata )[ which( som.nodes %in% nodes ) ]
			geneset.ids = unique( na.omit( gene.ids[ geneset.genes ] ) )
	
	
			GS.infos.underexpression[[ LETTERS[i] ]] = list()
	
			GS.infos.underexpression[[ LETTERS[i] ]]$metagenes = as.numeric( nodes )
			GS.infos.underexpression[[ LETTERS[i] ]]$genes = geneset.genes
	
			GS.infos.underexpression[[ LETTERS[i] ]]$mask = rep( NA, (preferences$dim.som1*preferences$dim.som1) )
			GS.infos.underexpression[[ LETTERS[i] ]]$mask[ as.numeric( nodes ) ] = 1
	
			if( preferences$geneset.analysis )
				GS.infos.underexpression[[ LETTERS[i] ]]$HG.p = GeneSet.HG( geneset.ids, unique.protein.ids, gs.def.list, sort=T )$p.value
		}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
		##### Positive Peaks ######
	
	
	
		
		peaks = apply( metadata, 1, max )
	
	
	
		GS.infos.positivepeaks = list()
	
		GS.infos.positivepeaks$overview.map = peaks
	
	
	
		
		
	
	
		peaks = matrix( peaks, preferences$dim.som1, preferences$dim.som1 )
		peaks[ which( peaks < quantile( peaks, preferences$summary.spot.cutoff ) ) ] = NA
	
		e.cluster = matrix( NA, preferences$dim.som1, preferences$dim.som1 )
		e.cluster[ which( !is.na(peaks) ) ] = -1
		count.cluster = 1	
	
	
		while( any( !is.na( peaks ) ) )
		{
			
			start.pix = which( peaks == max(peaks,na.rm=T), arr.ind=T )[1,]
	
			e.cluster = col.pix( e.cluster, start.pix[1], start.pix[2], count.cluster )
	
			peaks[ which( e.cluster == count.cluster ) ] = NA
	
			count.cluster = count.cluster + 1
		}
	
	
	
		e.cluster = as.vector( e.cluster )
	
		for( i in 1:length( na.omit( unique( e.cluster ) ) ) )
		{
			nodes = which( e.cluster == i )
	
	
			geneset.genes = rownames( indata )[ which( som.nodes %in% nodes ) ]
			geneset.ids = unique( na.omit( gene.ids[ geneset.genes ] ) )
	
	
			GS.infos.positivepeaks[[ LETTERS[i] ]] = list()
	
			GS.infos.positivepeaks[[ LETTERS[i] ]]$metagenes = as.numeric( nodes )
			GS.infos.positivepeaks[[ LETTERS[i] ]]$genes = geneset.genes
	
			GS.infos.positivepeaks[[ LETTERS[i] ]]$mask = rep( NA, (preferences$dim.som1*preferences$dim.som1) )
			GS.infos.positivepeaks[[ LETTERS[i] ]]$mask[ as.numeric( nodes ) ] = 1
	
			if( preferences$geneset.analysis )
				GS.infos.positivepeaks[[ LETTERS[i] ]]$HG.p = GeneSet.HG( geneset.ids, unique.protein.ids, gs.def.list, sort=T )$p.value
		}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
		##### Negative Peaks ######
	
	
	
		
		peaks = apply( metadata, 1, min )
	
	
	
		GS.infos.negativepeaks = list()
	
		GS.infos.negativepeaks$overview.map = peaks
	
	
	
		
		
	
	
		peaks = matrix( peaks, preferences$dim.som1, preferences$dim.som1 )
		peaks[ which( peaks > quantile( peaks, 1-preferences$summary.spot.cutoff ) ) ] = NA
	
		e.cluster = matrix( NA, preferences$dim.som1, preferences$dim.som1 )
		e.cluster[ which( !is.na(peaks) ) ] = -1
		count.cluster = 1	
	
	
		while( any( !is.na( peaks ) ) )
		{
			
			start.pix = which( peaks == min(peaks,na.rm=T), arr.ind=T )[1,]
	
			e.cluster = col.pix( e.cluster, start.pix[1], start.pix[2], count.cluster )
	
			peaks[ which( e.cluster == count.cluster ) ] = NA
	
			count.cluster = count.cluster + 1
		}
	
	
	
		e.cluster = as.vector( e.cluster )
	
		for( i in 1:length( na.omit( unique( e.cluster ) ) ) )
		{
			nodes = which( e.cluster == i )
	
	
			geneset.genes = rownames( indata )[ which( som.nodes %in% nodes ) ]
			geneset.ids = unique( na.omit( gene.ids[ geneset.genes ] ) )
	
	
			GS.infos.negativepeaks[[ LETTERS[i] ]] = list()
	
			GS.infos.negativepeaks[[ LETTERS[i] ]]$metagenes = as.numeric( nodes )
			GS.infos.negativepeaks[[ LETTERS[i] ]]$genes = geneset.genes
	
			GS.infos.negativepeaks[[ LETTERS[i] ]]$mask = rep( NA, (preferences$dim.som1*preferences$dim.som1) )
			GS.infos.negativepeaks[[ LETTERS[i] ]]$mask[ as.numeric( nodes ) ] = 1
	
			if( preferences$geneset.analysis )
				GS.infos.negativepeaks[[ LETTERS[i] ]]$HG.p = GeneSet.HG( geneset.ids, unique.protein.ids, gs.def.list, sort=T )$p.value
		}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
		##### Correlation Cluster ######
	
	
		GS.infos.correlation = list()
	
	
	
		c.map = cor( t( metadata ) )
		diag( c.map ) = NA
		rownames( c.map ) = c( 1:( preferences$dim.som1*preferences$dim.som1 ) )
		colnames( c.map ) = c( 1:( preferences$dim.som1*preferences$dim.som1 ) )
	
		c.cluster = rep( NA, (preferences$dim.som1*preferences$dim.som1) )
		names( c.cluster ) = c( 1:( preferences$dim.som1*preferences$dim.som1 ) )
	
	
	
		count.cluster = 1
	
	
		while( is.matrix(c.map) && nrow(c.map) > 0 && ncol(c.map) > 0  )
		{
	
			start.node = rownames(c.map)[ which( c.map == max( c.map, na.rm=T ), arr.ind=T )[1,1] ]
	
			cluster = names( which( c.map[ start.node, ] > 0.9 ) )
			cluster = c( start.node, cluster )
		
			if( length(cluster) > 1 )
			{
				c.cluster[ cluster ] = count.cluster
	
				geneset.genes = rownames( indata )[ which( som.nodes %in% as.numeric( cluster ) ) ]
				geneset.ids = unique( na.omit( gene.ids[ geneset.genes ] ) )
	
	
				GS.infos.correlation[[ count.cluster ]] = list()
	
				GS.infos.correlation[[ count.cluster ]]$metagenes = as.numeric( cluster)
				GS.infos.correlation[[ count.cluster ]]$genes = geneset.genes
	
				GS.infos.correlation[[ count.cluster ]]$mask = rep( NA, (preferences$dim.som1*preferences$dim.som1) )
				GS.infos.correlation[[ count.cluster ]]$mask[ as.numeric( cluster) ] = 1
	                           
				if( preferences$geneset.analysis )
					GS.infos.correlation[[ count.cluster ]]$HG.p = GeneSet.HG( geneset.ids, unique.protein.ids, gs.def.list, sort=T )$p.value
	
			
				count.cluster = count.cluster + 1
			}
	
			c.map = c.map[ -which( rownames(c.map)%in%cluster), -which( colnames(c.map)%in%cluster) ]
	
	
		}
		GS.infos.correlation$overview.map = c.cluster
	
	
	
	
	
	
	
	
	
	
	
		##### K-Means Clustering #####
	
	
		
		v = 0
		n.cluster = 1
		while( n.cluster < 100 && v < 0.8 )			# increase n.clusters til explained variance >= 80%
		{
			n.cluster = n.cluster + 1
	
			res = kmeans( metadata, n.cluster )
	
			if( exists( "res$betweenss" ) && exists(res$totss) )		# catch old R version (kmeans provides less values)
	
				v=res$betweenss / res$totss 
	
			else
			{
				v = 1
				n.cluster = 10
			}
		}
	
	
	
		GS.infos.kmeans = list()
	
		GS.infos.kmeans$overview.map = matrix( res$cluster, preferences$dim.som1, preferences$dim.som1 )
	
	
		for( i in 1:n.cluster )
		{
			nodes = which( res$cluster == i )
	
	
			geneset.genes = rownames( indata )[ which( som.nodes %in% nodes ) ]
			geneset.ids = unique( na.omit( gene.ids[ geneset.genes ] ) )
	
	
			GS.infos.kmeans[[ LETTERS[i] ]] = list()
	
			GS.infos.kmeans[[ LETTERS[i] ]]$metagenes = as.numeric( nodes )
			GS.infos.kmeans[[ LETTERS[i] ]]$genes = geneset.genes
	
			GS.infos.kmeans[[ LETTERS[i] ]]$mask = rep( NA, (preferences$dim.som1*preferences$dim.som1) )
			GS.infos.kmeans[[ LETTERS[i] ]]$mask[ as.numeric( nodes ) ] = 1
	
			if( preferences$geneset.analysis )
				GS.infos.kmeans[[ LETTERS[i] ]]$HG.p = GeneSet.HG( geneset.ids, unique.protein.ids, gs.def.list, sort=T )$p.value
		}
	




	cat( "Plotting Sample Profiles\n\n" ); flush.console()
	
	### Profile-Catalog
	
	
		pdf( paste( files.name, "- Results/Expression Profiles.pdf" ) , 29.7/2.54, 21/2.54 )
	
	
		par( mfrow=c( 7, 12 ) )
		par( mar=c(0.3,0.9,4.5,0.9) )
		count.col = 0
		for( gl in 1:length( unique( group.labels ) ) )
		{
			plot( 0, type="n", axes=F, xlab="", ylab="", xlim=c(0,1) )
			if( length( unique( group.labels ) ) > 1 )
				mtext( unique( group.labels )[gl], side=3, line = 2, cex=1.5, at=0, font=3, adj=0, col=unique(group.colors)[gl] )
	
			par(new=T)
			for( j in which( group.labels == unique( group.labels )[gl] ) )
			{
				image( matrix( metadata[,j], preferences$dim.som1, preferences$dim.som1 ), axes=F, col = colramp(1000) )
				title( paste(j,":",colnames(indata)[j]), line=1, cex.main=0.8 )
				box()
				
	
				count.col = count.col + 1 
			}
			if( count.col %% 12 != 0 )
			for( j in 1:( 12 - count.col %% 12 ) )
			{
				plot( 0, type="n", axes=F, xlab="", ylab="" )
				
				count.col = count.col + 1
			}
		}
	
	
		dev.off()
	
	
	
	
	
	
	
		pdf( paste( files.name, "- Results/Expression Profiles alternative.pdf" ) , 29.7/2.54, 21/2.54 )
	
	
	
	
		par( mar=c(0,0,0,0), mfrow=c(1,1) )
		plot( 0, type="n", xlab="", ylab="", axes=T )
		text( 1, 0.3, "WAD Profiles", cex=2 )
	
	
	
		par( mfrow=c( 7, 12 ) )
		par( mar=c(0.3,0.9,4.5,0.9) )
		WADv = c()
		count.col = 0
		for( gl in 1:length( unique( group.labels ) ) )
		{
			plot( 0, type="n", axes=F, xlab="", ylab="", xlim=c(0,1) )
			if( length( unique( group.labels ) ) > 1 )
				mtext( unique( group.labels )[gl], side=3, line = 2, cex=1.5, at=0, font=3, adj=0, col=unique(group.colors)[gl] )
	
			par(new=T)
			for( j in which( group.labels == unique( group.labels )[gl] ) )
			{
	
	
				meta = metadata[,j]
				
				w = ( meta - min( meta ) ) / ( max( meta ) - min( meta ) )
	
	
				WAD = w * meta			
	
				image( matrix( WAD, preferences$dim.som1, preferences$dim.som1 ), axes=F, col = colramp(1000), cex.main=0.6 )
				title( paste(j,":",colnames(indata)[j]), line=1, cex.main=0.8 )
				box()
	
				WADv = c( WADv, WAD )
	
				count.col = count.col + 1 
			}
			if( count.col %% 12 != 0 )
			for( j in 1:( 12 - count.col %% 12 ) )
			{
				plot( 0, type="n", axes=F, xlab="", ylab="" )
				
				count.col = count.col + 1
			}
		}
	
	
	
	
	
	
	
	
	
	
		par( mar=c(0,0,0,0), mfrow=c(1,1) )
		plot( 0, type="n", xlab="", ylab="", axes=T )
		text( 1, 0.3, "loglog FC Profiles", cex=2 )
	
	
	
		par( mfrow=c( 7, 12 ) )
		par( mar=c(0.3,0.9,4.5,0.9) )
		count.col = 0
		for( gl in 1:length( unique( group.labels ) ) )
		{
			plot( 0, type="n", axes=F, xlab="", ylab="", xlim=c(0,1) )
			if( length( unique( group.labels ) ) > 1 )
				mtext( unique( group.labels )[gl], side=3, line = 2, cex=1.5, at=0, font=3, adj=0, col=unique(group.colors)[gl] )
	
			par(new=T)
			for( j in which( group.labels == unique( group.labels )[gl] ) )
			{		
				meta = metadata[,j]
				meta.sign = sign( meta )
				meta = log10( abs( meta ) )
				meta = meta - min( meta, na.rm=T )
				meta = meta * meta.sign
	
				image( matrix( meta, preferences$dim.som1, preferences$dim.som1 ), axes=F, col = colramp(1000) )
				title( paste(j,":",colnames(indata)[j]), line=1, cex.main=0.8 )
				box()
				
	
				count.col = count.col + 1 
			}
			if( count.col %% 12 != 0 )
			for( j in 1:( 12 - count.col %% 12 ) )
			{
				plot( 0, type="n", axes=F, xlab="", ylab="" )
				
				count.col = count.col + 1
			}
		}
	
	
	
	
		dev.off()




	cat( "2nd level SOM\n\n" ); flush.console()


	supersom.custom = som( t(metadata), xdim=preferences$dim.som2, ydim=preferences$dim.som2 )
	
	
	
		if( preferences$dim.som2 == 20 )
		{
			supersom.20 = supersom.custom
		} else
		{
			supersom.20 = som( t(metadata), xdim=20, ydim=20 )
		}
	
	
	
	
	
		pdf( paste( files.name, "- Results/2nd level SOM.pdf" ), 21/2.54, 21/2.54 )
	
	
	
		##### Plot Supersom #####
	
		par( mar=c(1,1,1,1) )
	
	
	
	
		xl = c( min( supersom.custom$visual[,"x"] )-0.2, max( supersom.custom$visual[,"x"] )+0.2 )
		yl = c( min( -supersom.custom$visual[,"y"] )-0.2, max( -supersom.custom$visual[,"y"] )+0.2 )
		plot( supersom.custom$visual[,"x"], -supersom.custom$visual[,"y"], type="n", axes=F, xlab="", ylab="", xlim=xl, ylim=yl )
	
		if( ncol(indata) < 100 )
		{
			legend( "bottomright",  paste( 1:length( colnames( indata ) ), ":", colnames( indata ) ), cex=0.5, text.col=group.colors, ncol=(ncol(indata)-1)%/%25+1, bg="white" )
		}
		if( length(unique( group.labels ) ) > 1 )
		{
			legend( "topright", as.character(unique( group.labels )), cex=0.5, text.col=unique(group.colors), bg="white" )	
		}
	
		for( j in 1:nrow(supersom.custom$code.sum) )
		{
	
			which.samples = intersect( which( supersom.custom$visual[,"x"] == supersom.custom$code.sum[j,"x"] ),
							which( supersom.custom$visual[,"y"] == supersom.custom$code.sum[j,"y"] )	)
	
			if( !is.na( which.samples[1] ) )
			{
	
				which.samples = which.samples[ 1:min( 4, length(which.samples) ) ]
	
	
				x.seq = rep( c( -0.2, 0.2 ), 2 ) [ 1:length(which.samples) ]
				y.seq = c( rep( 0.2, 2 ), rep( -0.2, 2 ) ) [ 1:length(which.samples) ]
	
	
				points( supersom.custom$visual[ which.samples[1], "x" ]+x.seq, -supersom.custom$visual[ which.samples[1], "y" ]+y.seq,
					pch=16, col=group.colors[ which.samples ], cex=2.5 )
	
				points( supersom.custom$visual[ which.samples[1], "x" ]+x.seq, -supersom.custom$visual[ which.samples[1], "y" ]+y.seq,
					pch=1, col="black", cex=2.5, lwd=2 )
	
				text( supersom.custom$visual[ which.samples[1], "x" ]+x.seq, -supersom.custom$visual[ which.samples[1], "y" ]+y.seq,
					which.samples, col="black", cex=0.9 )
							
			}
		}
	
		box()
	
	
	
	
	
		##### Plot Supersom with real expression profiles ######	
	
	
	
		par( mar=c(1,1,1,1) )
	
	
	
		xl = c( min( supersom.20$visual[,"x"] )-0.2, max( supersom.20$visual[,"x"] )+0.2 )
		yl = c( min( -supersom.20$visual[,"y"] )-0.2, max( -supersom.20$visual[,"y"] )+0.2 )
		plot( supersom.20$visual[,"x"], -supersom.20$visual[,"y"], type="p", axes=F, xlab="", ylab="", xlim=xl, ylim=yl )
	
		if( ncol(indata) < 100 )
		{
			legend( "bottomright",  paste( 1:length( colnames( indata ) ), ":", colnames( indata ) ), cex=0.5, text.col=group.colors, ncol=(ncol(indata)-1)%/%25+1, bg="white" )
		}
		if( length(unique( group.labels ) ) > 1 )
		{
			legend( "topright", as.character(unique( group.labels )), cex=0.5, text.col=unique(group.colors), bg="white" )	
		}
	
		for( j in 1:nrow(supersom.20$code.sum) )
		{
	
			which.samples = intersect( which( supersom.20$visual[,"x"] == supersom.20$code.sum[j,"x"] ),
							which( supersom.20$visual[,"y"] == supersom.20$code.sum[j,"y"] )	)
	
			if( !is.na( which.samples[1] ) )
			{
	
				m = matrix( metadata[, which.samples[1] ], preferences$dim.som1, preferences$dim.som1 )
				if( max(m) - min(m) != 0 ) 
					m = ( m - min(m) ) / ( max(m) - min(m) )*999
				m = cbind( apply( m, 1, function(x){x} ) )[ nrow(m):1, ]
	
				x <- pixmapIndexed( m , col = colramp(1000) )
				addlogo( x, supersom.20$visual[ which.samples[1], "x" ]+c(-0.45,0.455), -supersom.20$visual[ which.samples[1], "y" ]+c(-0.45,0.45) )
	
	
				which.samples = which.samples[ 1:min( 4, length(which.samples) ) ]
	
	
				x.seq = rep( c( -0.2, 0.2 ), 2 ) [ 1:length(which.samples) ]
				y.seq = c( rep( 0.2, 2 ), rep( -0.2, 2 ) ) [ 1:length(which.samples) ]
	
	
				points( supersom.20$visual[ which.samples[1], "x" ]+x.seq, -supersom.20$visual[ which.samples[1], "y" ]+y.seq,
					pch=16, col=group.colors[ which.samples ], cex=2.5 )
	
				points( supersom.20$visual[ which.samples[1], "x" ]+x.seq, -supersom.20$visual[ which.samples[1], "y" ]+y.seq,
					pch=1, col="black", cex=2.5, lwd=2 )
	
				text( supersom.20$visual[ which.samples[1], "x" ]+x.seq, -supersom.20$visual[ which.samples[1], "y" ]+y.seq,
					which.samples, col="black", cex=0.9 )
							
			}
		}
	  
		box()
	
	
	
		dev.off()




	cat( "Processing Supporting Information\n\n" ); flush.console()
	

	pdf( paste( files.name, "- Results/Supporting Maps.pdf" ), 21/2.54, 21/2.54 )
	
	
	
	
	     	##### Plot Supporting Maps ######
	
	
		### Population Map ###
	
		par( mfrow = c(1,1) );	par( mar=c(5,6,4,5) )
		m = log10(som.result$code.sum[,"nobs"])
		m[ which( is.infinite(m) ) ] = 0
		image( matrix( m, preferences$dim.som1, preferences$dim.som1 ), axes=F, col = colramp(1000), main="Population Map", cex.main=2.5 )
			mtext( "log ( # genes in metagene )", side=1, line=1, cex=1.4 )
			box()
	
		par(new=T, mar=c(1,0,0,0) );	layout( matrix( c(0,0,0,0,1,0,0,0,0), 3, 3 ), c( 1, 0.05, 0.02 ), c( 0.14, 0.3, 1 ) )
		image( matrix( 1:100, 1, 100 ), col = colramp(1000), axes=F )
			axis( 2, at=c(0,1), c( min(som.result$code.sum[,"nobs"]), max(som.result$code.sum[,"nobs"]) ), las=2, tick=F, pos=-0.5, cex.axis=1.4 )	
			box()
		
	
	
	
	
	
		### Variance Map ###
	
		par( mfrow = c(1,1) );	par( mar=c(5,6,4,5) )
		image( matrix( log10( apply( metadata, 1, var ) ), preferences$dim.som1, preferences$dim.som1 ), axes=F, col = colramp(1000), main="Metagene Variance Map", cex.main=2.5 )
			mtext( "log ( metagene variance )", side=1, line=1, cex=1.4 )
			box()
		
		par(new=T, mar=c(1,0,0,0) );	layout( matrix( c(0,0,0,0,1,0,0,0,0), 3, 3 ), c( 1, 0.05, 0.02 ), c( 0.14, 0.3, 1 ) )
		image( matrix( 1:100, 1, 100 ), col = colramp(1000), axes=F )
			axis( 2, at=c(0,1), c( round(min(apply( metadata, 1, var ),na.rm=T),2), round(max(apply( metadata, 1, var ),na.rm=T),2) ), las=2, tick=F, pos=-0.5, cex.axis=1.4 )		
			box()
	
	
	
	
	
		### Covariance Map ###
	
		errors = c()
		for( i in 1:preferences$dim.som1^2 )
		{	
			genes = names( which( som.nodes == i ) )
	
			mean.cor = 0
			for( ii in genes )
			{
				suppressWarnings( {   	mean.cor = mean.cor + cor( metadata[ i, ], as.numeric(indata[ ii, ]) )   	} )
			}
			errors[i] = mean.cor / length( genes )
		}
	
		par( mfrow = c(1,1) );	par( mar=c(5,6,4,5) )
		image( matrix( errors, preferences$dim.som1, preferences$dim.som1 ), axes=F, col = colramp(1000), main="Gene-Metagene Covariance Map", cex.main=2.5 )
			mtext( "correlation genes - metagene", side=1, line=1, cex=1.4 )
			box()
	
		par(new=T, mar=c(1,0,0,0) );	layout( matrix( c(0,0,0,0,1,0,0,0,0), 3, 3 ), c( 1, 0.05, 0.02 ), c( 0.14, 0.3, 1 ) )
		image( matrix( 1:100, 1, 100 ), col = colramp(1000), axes=F )
			axis( 2, at=c(0,1), c( round(min(errors,na.rm=T),2),round(max(errors,na.rm=T),2) ), las=2, tick=F, pos=-0.5, cex.axis=1.4 )		
			box()
	
	
	
	
		### Deviation Map ###
	
		errors = rep( NA, preferences$dim.som1^2 )
		for( i in 1:preferences$dim.som1^2 )
		{	
			genes = names( which( som.nodes == i ) )
			e = NA
	
			if( length( genes ) > 1 )
			{
				e = apply( indata[ genes, ], 1, function(x){ 1/(ncol(indata)-1) * sum( ( x - metadata[ i, ] )^2 ) } )
			} else
			if( length( genes ) == 1 )
			{
				e = 1/(ncol(indata)-1) * sum( ( indata[ genes, ] - metadata[ i, ] )^2 )
			}	
		
			errors[i] = sqrt( mean(e) )
	
		}
	
		par( mfrow = c(1,1) );	par( mar=c(5,6,4,5) )
		image( matrix( errors, preferences$dim.som1, preferences$dim.som1 ), axes=F, col = colramp(1000), main="Gene-Metagene Deviation Map", cex.main=2.5 )
			mtext( "deviation genes - metagene", side=1, line=1, cex=1.4 )
			box()
	
	
		par(new=T, mar=c(1,0,0,0) );	layout( matrix( c(0,0,0,0,1,0,0,0,0), 3, 3 ), c( 1, 0.05, 0.02 ), c( 0.14, 0.3, 1 ) )
		image( matrix( 1:100, 1, 100 ), col = colramp(1000), axes=F )
			axis( 2, at=c(0,1), c( round(min(errors,na.rm=T),2), round(max(errors,na.rm=T),2) ), las=2, tick=F, pos=-0.5, cex.axis=1.4 )		
			box()
	
	
	
	
	
	
	
	
	
	     	##### Plot SOM 10x10 grid ######
	
		co.so = som.init( indata, xdim=10, ydim=10 )
	
	
		# Train SOM
	
		co.dim=10
		co.so = som.train( indata, co.so, xdim=10, ydim=10, 
					alpha=0.05, radius=10, rlen=nrow(indata)*2, inv.alp.c=nrow(indata)*2/100 )
		co.so = som.train( indata, co.so$code, xdim=10, ydim=10, 
					alpha=0.02, radius=min(3, 10), rlen=nrow(indata)*10, inv.alp.c=nrow(indata)*10/100 )
	
		co.so.code = co.so$code
		
	
		  
	
		par( mfrow=c(co.dim,co.dim) )
		par( mar=c(0,0,0,0) )
		for( j in (co.dim-1):0 )
		for( i in 0:(co.dim-1) )
		{
			ij = (i+1)+j*co.dim
	
	
			if( length( unique( group.labels) ) > 1 )
			{
				par( mar=c(5.4,0,0,0) )
				image( cbind(1:length(group.colors)), col=group.colors, axes=F ) 
				par(new=T, mar=c(0,0,0,0) )
			}
	
			plot( co.so.code[ij,], ylim=c(min(co.so.code),max(co.so.code)*1.2), xlab="", ylab="", axes=F, type="l", lwd=1.4 )		
			points( co.so.code[ij,], pch=16, cex=1 )
	
			abline( h=0, col="darkgray", lty=2, lwd=2 )
	
			text( 0.5, max(co.so.code)*0.65, co.so$code.sum[ ij, "nobs" ], cex=1.2, adj=0 )
	
			box()
	
		}
	
	
	
		dev.off()



	cat( "Processing 2nd level Standard Analyses\n\n" ); flush.console()
	


	
		## Heatmap Wrapper
		heatmap.wrap = function( RowSideColors, ColSideColors, ... )
		{
			if( !missing( RowSideColors ) && !missing( ColSideColors ) )
			{
				if( length( unique( RowSideColors ) ) > 1 && length( unique( ColSideColors ) ) > 1 )
				{
					heatmap( RowSideColors=RowSideColors, ColSideColors=ColSideColors, ... )
					return( NULL )
				}
			}else
			if( !missing( ColSideColors ) )
			{
				if( length( unique( ColSideColors ) ) > 1 )
				{
					heatmap( ColSideColors=ColSideColors, ... )
					return( NULL )
				}
			}	
			
			heatmap( ... )
		}
	
	
	
	
		filter.list = list()
	
		filter.list[[1]] = list( s = c(1:preferences$dim.som1^2), n = paste("all", preferences$dim.som1^2, "Metagenes") )
	
	
		if( preferences$dim.som1^2 > 1000 ) 
			filter.list[[ length(filter.list) + 1 ]] = list( s = order( apply( abs(metadata), 1, max ), decreasing=T )[1:1000], n = "High Expression: 1000 Metagenes" )
		if( preferences$dim.som1^2 > 100 ) 
			filter.list[[ length(filter.list) + 1 ]] = list( s = order( apply( abs(metadata), 1, max ), decreasing=T )[1:100], n = "High Expression: 100 Metagenes" )
	
	
		if( preferences$dim.som1^2 > 1000 ) 
			filter.list[[ length(filter.list) + 1 ]] = list( s = order( apply( abs(metadata), 1, var ), decreasing=T )[1:1000], n = "Variance: 1000 Metagenes" )
		if( preferences$dim.som1^2 > 100 ) 
			filter.list[[ length(filter.list) + 1 ]] = list( s = order( apply( abs(metadata), 1, var ), decreasing=T )[1:100], n = "Variance: 100 Metagenes" )
	
	
	
	
	
	
	
	
	
	
		pdf( paste( files.name, "- Results/2nd level correlation analysis.pdf" ), 29.7/2.54, 21/2.54 )
	
	
		for( i in 1:length(filter.list) )
		{
	
			s = metadata[ filter.list[[i]]$s , ]
			par( mar=c(1,1,1,1) )
	
	
	
			adj.matrix = cor( metadata ) * -1
			g = graph.adjacency( adj.matrix, weighted=T,  mode="undirected" )
			stg = minimum.spanning.tree( g ) 
	
			layout=layout.kamada.kawai(stg)
	
			plot( stg, layout=layout, vertex.size=5, vertex.label = colnames(indata), vertex.label.cex=1.4, vertex.color=group.colors, main=filter.list[[i]]$n )
				box()
	
			plot( stg, layout=layout, vertex.size=5, vertex.label = rep("",ncol(indata)) , vertex.label.cex=1.4, vertex.color=group.colors, main=filter.list[[i]]$n )
				box()
	
	
	
	
	
	
			heatmap.wrap( x=cor( s ), col=colramp(1000), scale="n", main=filter.list[[i]]$n, mar=c(8,8), ColSideColors=group.colors, RowSideColors=group.colors )
				par(new=T)
				plot(0,type="n", axes=F, xlab="", ylab="" )
	
			heatmap.wrap( x=cor( s ), col=colramp(1000), scale="n", main=filter.list[[i]]$n, mar=c(8,8), ColSideColors=group.colors, RowSideColors=group.colors, Colv=NA, Rowv=NA )
				par(new=T)
				plot(0,type="n", axes=F, xlab="", ylab="" )
	
	
	
		}
	
	
		dev.off()
	
	
	
	
	
	
	
	
	
	
	
		pdf( paste( files.name, "- Results/2nd level distance analysis.pdf" ), 29.7/2.54, 21/2.54 )
	
	
		for( i in 1:length(filter.list) )
		{
	
			s = metadata[ filter.list[[i]]$s , ]
			par( mar=c(1,1,1,1) ) 
	
	
	
			heatmap.wrap( x=s, col=colramp(1000), main=filter.list[[i]]$n, margins = c(10, 5), scale="n", labRow=NA, ColSideColors=group.colors )
				par(new=T)
				plot(0,type="n", axes=F, xlab="", ylab="" )
	
			heatmap.wrap( x=s, col=colramp(1000), main=filter.list[[i]]$n, margins = c(10, 5), scale="n", labRow=NA, ColSideColors=group.colors, Colv=NA )
				par(new=T)
				plot(0,type="n", axes=F, xlab="", ylab="" )
	
	
	
	
			if( ncol(metadata) > 2 )
			{	
				phylo.tree = nj( dist( t( s ) ) )
				phylo.tree$tip.label = colnames(indata)
			
				plot.phylo( phylo.tree, "unrooted", cex=0.5, tip.color=group.colors )
					title(filter.list[[i]]$n)
					box()
			}
	
	
		}
	
	
		dev.off()
	
	
	
	
	
	
	
	
	
	
		pdf( paste( files.name, "- Results/2nd level component analysis.pdf" ), 29.7/2.54, 21/2.54 )
	
	
		for( i in 1:length(filter.list) )
		try(
		{
	
			s = metadata[ filter.list[[i]]$s , ];
			par( mar=c(1,1,1,1) );
	
	
	
			ICA.metagenes = fastICA( t( s ), 2 )$S;
	
	
			plot( ICA.metagenes[,1], ICA.metagenes[,2], type="p", pch=16, col=group.colors, cex=4, axes=F, xlab="", ylab="", main=filter.list[[i]]$n );
				if( ncol(indata) < 100 )			{ legend( "bottomright",  paste( 1:length( colnames( indata ) ), ":", colnames( indata ) ), cex=0.5, text.col=group.colors, ncol=(ncol(indata)-1)%/%25+1, bg="white" )	};
				if( length(unique( group.labels ) ) > 1 )	{ legend( "topright", as.character(unique( group.labels )), cex=0.5, text.col=unique(group.colors), bg="white" )	};
	
				points( ICA.metagenes[,1], ICA.metagenes[,2], pch=16, col=group.colors, cex=4 );
				points( ICA.metagenes[,1], ICA.metagenes[,2], pch=1, col="black", cex=4 );
				text( ICA.metagenes[,1], ICA.metagenes[,2], colnames(indata), col="white" );   #1:nrow(ICA.metagenes)
				box()
	
	
			plot( ICA.metagenes[,1], ICA.metagenes[,2], type="p", pch=16, col=group.colors, cex=4, axes=F, xlab="", ylab="", main=filter.list[[i]]$n );
				if( ncol(indata) < 100 )			{ legend( "bottomright",  paste( 1:length( colnames( indata ) ), ":", colnames( indata ) ), cex=0.5, text.col=group.colors, ncol=(ncol(indata)-1)%/%25+1, bg="white" )	};
				if( length(unique( group.labels ) ) > 1 )	{ legend( "topright", as.character(unique( group.labels )), cex=0.5, text.col=unique(group.colors), bg="white" )	};
	
				points( ICA.metagenes[,1], ICA.metagenes[,2], pch=16, col=group.colors, cex=4 );
				points( ICA.metagenes[,1], ICA.metagenes[,2], pch=1, col="black", cex=4 );
				box()
	
	
	
		} , silent=T )
	
	
		dev.off()




	cat( "Summary Sheets\n\n" ); flush.console()



	
		dir.create( paste( files.name, "- Results/Summary Sheets - Overviews" ), showWarnings=F )
	
	
	
	
	
		
	
	
	
		n.corr.cluster = 10 
		n.sets = 40
		n.samples = min( 20, ncol(indata) )
		
	
	
	
	
	
	
	
		plot.set.list = function( set.list, main )
		{
	
			if( main != "Correlation Cluster" )
			{
	
				layout( matrix( c( 1, 2 ), 1, 2 ), c( 2, 1 ), 1 )
				par( mar=c(5, 4, 4, 1) )
	
				image( matrix( set.list$overview.map, preferences$dim.som1, preferences$dim.som1 ), axes=F, col = colramp(1000), main=main, cex.main=1.5 )
					box()
		
				par( mar=c(5, 1, 4, 2) )
				plot( 0, type="n", axes=F, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i" )
					box()
	
			}
	
	
	
			# Spot Altas
	
	
			count = 1
			overview.map2 = list()
			overview.map2$mask = rep( NA, preferences$dim.som1^2 )
			overview.map2$top.GS.infos.samples = list()
			overview.map2$label.positions = list()
	
			for( m in which( names(set.list) != "overview.map" ) )
			{
				overview.map2$mask[ which( !is.na( set.list[[ m ]]$mask ) ) ] = count
				overview.map2$top.GS.infos.samples[[ count ]] = names( head( set.list[[ m ]]$HG.p , 3 ) )
	
				mid.x = mean( range( som.result$code.sum[ which( !is.na( set.list[[ m ]]$mask ) ), 1 ] ) )
				mid.y = mean( range( som.result$code.sum[ which( !is.na( set.list[[ m ]]$mask ) ), 2 ] ) )
				overview.map2$label.positions[[ count ]] = c( mid.x+0.5, mid.y+0.5 )
	
				count = count+1
			}
	
	
	
	
			layout( matrix( c( 1, 2 ), 1, 2 ), c( 2, 1 ), 1 )
			par( mar=c(5, 4, 4, 1) )
	
			image( matrix( overview.map2$mask, preferences$dim.som1, preferences$dim.som1 ), axes=F, col = colramp(count-1 ), main=main, cex.main=1.5 )
				box()
	
			par( new=T )
			plot( 0, type="n", axes=T, xlab="", ylab="", xlim=c(0,preferences$dim.som1), ylim=c(0,preferences$dim.som1), xaxs="i", yaxs="i" )
		
			points( do.call( rbind, overview.map2$label.positions ), pch=16, cex=3, col="black" )
			points( do.call( rbind, overview.map2$label.positions ), pch=1, cex=3, col="white" )
			text( do.call( rbind, overview.map2$label.positions ), LETTERS[ 1:(count-1) ], col="white" )
	
	
			par( mar=c(5, 1, 4, 2) )
			plot( 0, type="n", axes=F, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i" )
			box()
	
			if( preferences$geneset.analysis )
			{
				leg.col = colramp( count-1 )
				leg.col = as.vector( sapply( leg.col, c, NA, NA ) )
				leg.num = LETTERS[ 1:(count-1) ]
				leg.num = as.vector( sapply( leg.num, c, NA, NA ) )
	
				legend( x=0.05, y=1, unlist(overview.map2$top.GS.infos.samples), cex=0.7, col = leg.col, pch=15, pt.cex=1.5, bty="n"  )  
				legend( x=-0.04, y=1, legend=leg.num, cex=0.7, bty="n" ) 
			}
	
	
	
	
	
	
	
	
	
			# Spot - Sample - Heatmap
	
	
			sample.spot.expression = matrix( NA, 0, ncol(indata) )
	
			for( m in which( names(set.list) != "overview.map" ) )
			{
				mean.FC = apply( metadata, 2, function(x){ max( x[ set.list[[ m ]]$metagenes ] ) } )
	
				sample.spot.expression = rbind( sample.spot.expression, mean.FC )
			}
			rownames( sample.spot.expression ) = LETTERS[ 1:nrow(sample.spot.expression) ]
	
	
	
	
			layout( matrix( c( 0,2,0,3,1,0,0,4,5 ), 3, 3), heights=c(0.8,6,2), widths=c(0.5,5,3) )
		
			par( mar=c(0,0,0,0) )
	
		
			sample.spot.expression.image = if( nrow(sample.spot.expression) > 1 ) -t(sample.spot.expression[nrow(sample.spot.expression):1,]) else as.matrix( -(sample.spot.expression[nrow(sample.spot.expression):1,]) )
			
			image( 1:ncol(indata), 1:nrow(sample.spot.expression), sample.spot.expression.image, axes=F, ylim=0.5+c( 0,nrow(sample.spot.expression)), yaxs="i", xlab="", ylab="" )
				box()
	  			axis(1, 1:ncol(indata), labels = colnames(indata), las = 2, line = -0.5, tick = 0, cex.axis=1.4 )
	
	
			plot( 0, type="n", xlab="", ylab="", axes=F, xlim=c(0,1), ylim=0.5+c( 0,nrow(sample.spot.expression)), yaxs="i" )	
				text( 0.7, nrow(sample.spot.expression):1, rownames(sample.spot.expression), adj=1, cex=1.8 )
			
	
	
			par( mar=c(1,0,2,0) )
			plot( 0, type="n", xlab="", ylab="", axes=F )
			
	
			par( mar=c(0,0,0,0) )
			plot( 0, type="n", xlab="", ylab="", axes=F, xlim=c(0,1), ylim=0.5+c( 0,nrow(sample.spot.expression)), yaxs="i" )
				pos = as.vector( sapply( c(1:nrow(sample.spot.expression)), function(x){ c( x-0.26, x, x+0.26 )} ) )
				text( 0.05, rev(pos) , unlist(overview.map2$top.GS.infos.samples), adj=0, cex=1 ) #0.6
	
	
	
			par( mar=c(5,2,4,2) )
			image( matrix( 1:100, 100, 1 ), col=colorRampPalette( c("white","yellow","red") )(1000), axes=F, xlab="" )
				axis( 1, at=c(0,1), round( c( min(-sample.spot.expression.image), max(-sample.spot.expression.image) ), 1 ), las=2, tick=F, pos=-0.8, cex.axis=1.4 )	
				mtext( expression(paste("<",Delta,"e", '' ^ meta, ">")), side=1, line=0.5 )
				box()
		
	
	
	
	
	
	
			# Individual spot sheets
	
			for( m in which( names(set.list) != "overview.map" ) )
			{
	
				sample.with.spot = apply( metadata, 2, function(x)
							{  
								sample.spot = which( x > quantile( x, 0.98 ) ) 
								length( intersect( set.list[[ m ]]$metagenes, sample.spot ) ) > 0
							} )
				
	
	
	
				layout( matrix( c(1,2,0,1,3,0,4,4,5,6,6,7), 3, 4 ), width=c(1,1,2,2), heights=c(2,1,1) )
	
	
	
	
				par( mar=c( 0,0,0,0 ) )
	
				plot( 0, type="n", axes=F, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1) )
	
				text( 0.1, 0.94, main , cex=2.6, adj=0 )	
	
		 		text( 0.1, 0.8, "Spot Summary" , cex=1.8, adj=0 )	
	
				text( 0.1, 0.7, paste( "# metagenes =", length( set.list[[ m ]]$metagenes ) ), adj=0 )
				text( 0.1, 0.65, paste( "# genes =", length( set.list[[ m ]]$genes ) ), adj=0 )
	
	
				text( 0.1, 0.55, paste( "<r> metagenes =", round(mean( cor( t( metadata[ set.list[[ m ]]$metagenes, ] ) ) ), 2 ) ), adj=0 )
				suppressWarnings( try( text( 0.1, 0.5, paste( "<r> genes =", round(mean( cor( t( indata[ set.list[[ m ]]$genes, ] ) ) ), 2 ) ), adj=0 ), silent=T ) )
	
				text( 0.1, 0.4, paste( "# samples with spot =", sum(sample.with.spot), "(", round( 100 * sum(sample.with.spot)/ncol(indata), 1), "% )" ), adj=0 )
	
	
	
	
	
	
	
				par( mar=c( 2,3,3,1 ) )
	
				image( matrix( set.list$overview.map, preferences$dim.som1, preferences$dim.som1 ), axes=F, col = colramp(1000), main="Overview Map", cex.main=1.5 )
					axis( 1, seq( 0, 1, length.out = preferences$dim.som1/10+1 ), c( 1, seq( 10, preferences$dim.som1, length.out = preferences$dim.som1/10 ) ), cex.axis=1.0 )
					axis( 2, seq( 0, 1, length.out = preferences$dim.som1/10+1 ), c( 1, seq( 10, preferences$dim.som1, length.out = preferences$dim.som1/10 ) ), cex.axis=1.0, las=1 )
					box()
	
				image( matrix( set.list$overview.map, preferences$dim.som1, preferences$dim.som1 ), axes=F, col = colramp(1000), main="Spot", cex.main=1.5 )
					par( new=T )
					mask = set.list[[ m ]]$mask
					mask[ which( is.na( set.list[[ m ]]$mask ) ) ] = 1
					mask[ which( !is.na( set.list[[ m ]]$mask ) ) ] = NA
					image( matrix( mask, preferences$dim.som1, preferences$dim.som1 ), axes=F, col = "white" )
					axis( 1, seq( 0, 1, length.out = preferences$dim.som1/10+1 ), c( 1, seq( 10, preferences$dim.som1, length.out = preferences$dim.som1/10 ) ), cex.axis=1.0 )
					axis( 2, seq( 0, 1, length.out = preferences$dim.som1/10+1 ), c( 1, seq( 10, preferences$dim.som1, length.out = preferences$dim.som1/10 ) ), cex.axis=1.0, las=1 )
					box()
	 
	
	
	
	
	
	
	
	
				par( mar=c( 0,0,0,0 ) )
	
				x.coords = c( 0, 0.15, 0.4 )
				y.coords = seq( 0.75, 0.02, length.out=n.samples )
	
	
			
	
	
	
	
				mean.FC = apply( indata, 2, function(x){ mean( x[ set.list[[ m ]]$genes ] ) } )
	
				o = order( abs( mean.FC ), decreasing=T )[1:n.samples]
	
				plot( 0, type="n", axes=F, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1) )
			
	
					text( 0, 0.88, "Sample List", cex=1.8, adj=0 )
	
	
					text( x.coords, rep( 0.82, 3 ), c( "Rank", "Sample", "<log(FC)>" ), cex=1, adj=0 )
			
					text( x.coords[1], y.coords, c( 1:n.samples ), adj=0 )
	
					text( x.coords[2], y.coords, colnames( indata )[o], cex=0.6, adj=0 )
					text( x.coords[3], y.coords, round( mean.FC[o], 2 ), cex=0.6, adj=0 )
			
	
	
	
				plot( 0, type="n", axes=F, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1) )
	
	
	
	
				if( preferences$geneset.analysis )
				{
		
					top.gs.p = sort( set.list[[ m ]]$HG.p )[ 1:n.sets ]
	
	
					par( mar=c( 0,0,0,0 ) )
	
					x.coords = c( 0, 0.1, 0.23, 0.31 )
					y.coords = seq( 0.75, 0.02, length.out=n.sets )
	
	
	
					plot( 0, type="n", axes=F, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1) )
	
						text( 0, 0.88, "Geneset Overrepresentation", cex=1.8, adj=0 )
	
						text( x.coords, 0.82, c( "Rank", "p-value", "Geneset", "" ), cex=1, adj=0 )
			
						text( x.coords[1], y.coords, c( 1:n.sets ), adj=0 )
						text( x.coords[2], y.coords, format( top.gs.p, digits=1 ), cex=0.6, adj=0 )
						text( x.coords[3], y.coords, sapply( gs.def.list, function(x){ x$Type } )[names(top.gs.p)], cex=0.6, adj=0 )
						text( x.coords[4], y.coords, names( top.gs.p ), cex=0.6, adj=0 )
	
	
	
				}
	
	
			}
		
		
		}
	
	
	
	
		pdf( paste( files.name, " - Results/Summary Sheets - Overviews/Overexpression.pdf", sep="" ), 29.7/2.54, 21/2.54 )
	
		plot.set.list( set.list=GS.infos.overexpression, main="Sample-Overexpression" )
	
		dev.off()
	
	
		pdf( paste( files.name, " - Results/Summary Sheets - Overviews/Underexpression.pdf", sep="" ), 29.7/2.54, 21/2.54 )
	
		plot.set.list( GS.infos.underexpression, main="Sample-Underexpression" )
	
		dev.off()
	
	
		pdf( paste( files.name, " - Results/Summary Sheets - Overviews/Positive Metagenepeaks.pdf", sep="" ), 29.7/2.54, 21/2.54 )
	
		plot.set.list( GS.infos.positivepeaks, main="Metagene Maxima" )
	
		dev.off()
	
	
		pdf( paste( files.name, " - Results/Summary Sheets - Overviews/Negative Metagenepeaks.pdf", sep="" ), 29.7/2.54, 21/2.54 )
	
		plot.set.list( GS.infos.negativepeaks, main="Metagene Minima" )
	
		dev.off()
	
	
		pdf( paste( files.name, " - Results/Summary Sheets - Overviews/Correlation Cluster.pdf", sep="" ), 29.7/2.54, 21/2.54 )
	
		plot.set.list( GS.infos.correlation, main="Correlation Cluster" )
	
		dev.off()
	





	cat( "HTML Report\n\n" ); flush.console()
	


	
	
	outfile = file( paste( files.name,".html",sep=""), "w" )
	          	
	cat( "
		<html> <head> <TITLE>Summary of",files.name,"dataset</TITLE>
		</head> <body bgcolor=#FFFFFF >
		<style type=text/css> p{ margin-top: 1px; margin-bottom: 1px; padding-left: 10px; text-indent: -10px }
		</style> ", file = outfile)
	
	
	
	
	cat( "<H1>General Information</H1>
		<TABLE BORDER=2, WIDTH=90%>
		<colgroup>
			<col width=\"50%\">
			<col width=\"50%\">
		</colgroup>
	
		<TR>
		<TD>Number of samples:</TD>
		<TD>", ncol(indata), "</TD>
		</TR>
	
		<TR>
		<TD>Number of categories:</TD>
		<TD>", length(unique(group.labels)), "</TD>
		</TR>
	
		<TR>
		<TD>Number of genes:</TD>
		<TD>", nrow(indata), "</TD>
		</TR>
	
		<TR>
		<TD>Analysis finished:</TD>
		<TD>", format(Sys.time(), "%a %b %d %X %Y %Z"), "</TD>
		</TR>
	
		</TABLE><br>", sep="", file = outfile)
	
	
	
	
	cat( "<TABLE BORDER=2, WIDTH=90%>
		<colgroup>
			<col width=\"50%\">
			<col width=\"50%\">
		</colgroup>
	
		<TR>
		<TD>Dimension 1st level SOM:</TD>
		<TD>", preferences$dim.som1, "x", preferences$dim.som1, "</TD>
		</TR>
	
		<TR>
		<TD>Dimension 2nd level SOM:</TD>
		<TD>", preferences$dim.som2, "x", preferences$dim.som2, "</TD>
		</TR>
	
		</TABLE>", sep="", file = outfile)
	
	
	
	
	
	
	
	cat( "<H1>Results</H1>
	
		<TABLE BORDER=2 , WIDTH=90%>
		<colgroup>
			<col width=\"50%\">
			<col width=\"50%\">
		</colgroup>
		<thead><tr>
			<th colspan=2, BGCOLOR=\"#99CCFF\" >Maps (experiment atlas)</th>
		</tr>	</thead>
	
		
		<TR>
		<TD rowspan=5 >These reports show the collection of first level SOM of all tissue samples, supporting maps and the second level SOM. First level SOM are shown with different contrast (log FC-, WAD- and double log-scale).</TD>
	
		<TD><a href=\"", files.name, " - Results/Expression Profiles.pdf\" target=\"_blank\">
			<b>First level SOM expression profiles (FC)</b></a></TD>
		</TR>
	
		<TR>
		<TD><a href=\"", files.name, " - Results/Expression Profiles alternative.pdf\" target=\"_blank\">
			First level SOM expression profiles (WAD & loglog)</a></TD>
		</TR>
	
		<TR>
		<TD><a href=\"", files.name, " - Results/Supporting Maps.pdf\" target=\"_blank\">
			Supporting Maps</a></TD>
		</TR>
	
		<TR>
		<TD><a href=\"", files.name, " - Results/2nd level SOM.pdf\" target=\"_blank\">
			Second level SOM</a></TD>
		</TR>
	
		</TABLE><br>
	
	
	
	
	
		<TABLE BORDER=2, WIDTH=90%>
		<colgroup>
			<col width=\"50%\">
			<col width=\"50%\">
		</colgroup>
		<thead><tr>
			<th colspan=4, BGCOLOR=\"#99CCFF\" >Metagene based analysis</th>
		</tr>	</thead>
	
		
		<TR>
		<TD rowspan=5 >Several agglomerative methods based either on distance or on correlation metrics are applied to the samples using filtered subsets of metagenes.
				The reports show two-way hierarchical clustering heatmaps, pairwise correlation maps, minimum spanning trees and the ICA results.<br></TD>
	
		<TD><a href=\"", files.name, " - Results/2nd level distance analysis.pdf\" target=\"_blank\">
			Distance based methods ( Clustering, Phylogenetic Tree )</a></TD>
		</TR>
	
		<TR>
		<TD><a href=\"", files.name, " - Results/2nd level correlation analysis.pdf\" target=\"_blank\">
			Correlation based methods ( MST, PCM )</a></TD>
		</TR>
	
		<TR>
		<TD><a href=\"", files.name, " - Results/2nd level component analysis.pdf\" target=\"_blank\">
			Component based methods ( 2d-ICA )</a></TD>
		</TR>
	
	
		</TABLE><br>
	
	
	
	
	
		<TABLE BORDER=2, WIDTH=90%>
		<colgroup>
			<col width=\"50%\">
			<col width=\"50%\">
		</colgroup>
		<thead><tr>
			<th colspan=5, BGCOLOR=\"#99CCFF\" >Spot summaries</th>
		</tr>	</thead>
	
		
		<TR>
		<TD rowspan=5 >These analyses apply different criteria of spot selection such as overexpression, underexpression, maximum and minimum of metagene expression and mutual correlations between the metagenes.
				GO-enrichment analysis provides the three leading genes in the respective HG-enrichment list of each of the spots considered. Spot-related heatmaps characterize the expression profiles of the selected features in the series of samples.
				Single spot summary sheets provide detailed information on each of the spots such as the ranked list of samples which overexpress this feature in decreasing order according to the mean expression of the spot and the ranked list of the HG-overrepresented gene sets.</TD>
	
		<TD><a href=\"", files.name, " - Results/Summary Sheets - Overviews/Overexpression.pdf\" target=\"_blank\">
			Overexpression</a></TD>
		</TR>
	
		<TR>
		<TD><a href=\"", files.name, " - Results/Summary Sheets - Overviews/Underexpression.pdf\" target=\"_blank\">
			Underexpression</a></TD>
		</TR>
	
		<TR>
		<TD><a href=\"", files.name, " - Results/Summary Sheets - Overviews/Positive Metagenepeaks.pdf\" target=\"_blank\">
			Metagene Maxima</a></TD>
		</TR>
	
		<TR>
		<TD><a href=\"", files.name, " - Results/Summary Sheets - Overviews/Negative Metagenepeaks.pdf\" target=\"_blank\">
			Metagene Minima</a></TD>
		</TR>
	
		<TR>
		<TD><a href=\"", files.name, " - Results/Summary Sheets - Overviews/Correlation Cluster.pdf\" target=\"_blank\">
			Mutual Correlation</a></TD>
		</TR>
	
	
		</TABLE><br>", sep="", file = outfile)
	
	
	
	
	cat("	</body> </html> ", sep="", file = outfile)
	close(outfile)
	
	
	
	
	 



	cat( "Finished:", format(Sys.time(), "%a %b %d %X\n\n" ) )
	flush.console()


	return(invisible(0))
}




