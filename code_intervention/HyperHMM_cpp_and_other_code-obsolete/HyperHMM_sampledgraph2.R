## From
## https://github.com/rdiaz02/HT_test_run/blob/igj-hyperhmm/K_1_example/example_run.R

library(Rcpp)
library(RcppArmadillo)
library(igraph)
library(ggraph)
library(stringr)

## The next file must either exist or by a symlink to the
## existing one
## It comes from https://github.com/StochasticBiology/hypercube-hmm
sourceCpp("hyperhmm-r.cpp")


plotHyperHMM.sampledgraph2 = function(my.post, max.samps = 1e4, thresh = 0.05,
                                      node.labels = TRUE, use.arc = FALSE, no.times = FALSE,
                                      small.times = FALSE, times.offset = c(0.1,-0.1),
                                      edge.label.size = 2, edge.label.angle = "across",
                                      edge.label.colour = "#000000", edge.check.overlap = TRUE,
                                      featurenames = FALSE, truncate = -1,
                                      node.label.size = 2, use.timediffs = TRUE) {
    if(featurenames == TRUE) {
        featurenames = my.post$featurenames
    } else {
        featurenames = c("")
    }
    edge.from = edge.to = edge.time = edge.change = c()
    bigL = my.post$L
    if(truncate == -1 | truncate > bigL) { truncate = bigL }
    nsamps = min(max.samps, nrow(my.post$routes))
    ltraj = strsplit(my.post$viz, " ")
    for(i in 1:nsamps) {
        state1 = ltraj[[i]][1]
        state2 = ltraj[[i]][2]
        edge.from = c(edge.from, state1)
        edge.to = c(edge.to, state2)
        chars1 <- unlist(strsplit(state1, split = ""))
        chars2 <- unlist(strsplit(state2, split = ""))

                                        # Find positions where characters differ
        locus <- which(chars1 != chars2)

        edge.change = c(edge.change, locus)
    }

    df = data.frame(From=edge.from, To=edge.to, Change=edge.change)
    dfu = unique(df[,1:3])
    if(length(featurenames) > 1) {
        dfu$Change = featurenames[dfu$Change+1]
    }
    dfu$Flux = dfu$MeanT = dfu$SDT = NA
    for(i in 1:nrow(dfu)) {
        this.set = which(df$From==dfu$From[i] & df$To==dfu$To[i])
        dfu$Flux[i] = length(this.set)
        dfu$label[i] = paste(c("+", dfu$Change[i]), collapse="")


    }
    dfu$Flux = dfu$Flux / (nsamps/bigL)
    dfu = dfu[dfu$Flux > thresh,]
    trans.g = graph_from_data_frame(dfu)
                                        #bs = unlist(lapply(as.numeric(V(trans.g)$name), DecToBin, len=bigL))
                                        #bs = unlist(lapply(as.numeric(as.vector(V(trans.g))), DecToBin, len=bigL))
    bs = V(trans.g)$name
    V(trans.g)$binname = bs
    layers = str_count(bs, "1")

    if(truncate > bigL/2) {
        this.plot=  ggraph(trans.g, layout="sugiyama", layers=layers)
    } else {
        this.plot=  ggraph(trans.g, layout="tree")
    }
    if(use.arc == TRUE) {
        this.plot= this.plot +
            geom_edge_arc(aes(edge_width=Flux, edge_alpha=Flux, label=label, angle=45),
                          label_size = edge.label.size, label_colour=edge.label.colour, color="#AAAAFF",
                          label_parse = TRUE, angle_calc = edge.label.angle, check_overlap = edge.check.overlap) +
            scale_edge_width(limits=c(0,NA)) + scale_edge_alpha(limits=c(0,NA)) +
            theme_graph(base_family="sans")
    } else {
        this.plot=  this.plot +
            geom_edge_link(aes(edge_width=Flux, edge_alpha=Flux, label=label, angle=45),
                           label_size = edge.label.size, label_colour=edge.label.colour, color="#AAAAFF",
                           label_parse = TRUE, angle_calc = edge.label.angle, check_overlap = edge.check.overlap) +
            scale_edge_width(limits=c(0,NA)) + scale_edge_alpha(limits=c(0,NA)) +
            theme_graph(base_family="sans")
    }

    if(node.labels == TRUE) {
        this.plot = this.plot + geom_node_text(aes(label=binname),size=node.label.size)
    }

    return(this.plot)
}

# do the inference -- takes a couple seconds each
## res = HyperHMM(d_2000, nboot = 0)
## plotHyperHMM.sampledgraph2(res)
