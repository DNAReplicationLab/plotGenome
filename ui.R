library(shiny)
library(ggplot2)
library(colourpicker)
fluidPage(theme="plotGenome.css",
	headerPanel("plotGenome: web-based tool to analyse sortSeq data.", windowTitle="plotGenome"),
		
sidebarLayout(
	sidebarPanel(
			conditionalPanel(condition="input.tabs=='About'",
				div(id="aboutSide",
					HTML("<h4 style='padding-bottom:10px;'><b>Selected publications</b></h4><ul style='text-align:justify;'>
							<li>
								Müller, C. A., & Nieduszynski, C. A. (2012). <i>Conservation of replication timing reveals global and local 
								regulation of replication origin activity.</i> Genome Research, 22(10), 1953–1962.&nbsp;&nbsp;&nbsp;
								<a href='http://doi.org/10.1101/gr.139477.112' target='_blank'>View</a>&nbsp;&nbsp;&nbsp;
								<a href='https://www.ncbi.nlm.nih.gov/pubmed/22767388' target='_blank'>Pubmed</a>
							</li>
							<br>
							<li>
								Hawkins, M., Malla, S., Blythe, M. J., Nieduszynski, C. A., & Allers, T. (2013). 
								<i>Accelerated growth in the absence of DNA replication origins.</i> Nature, 503(7477), 544–547. 
								<a href='http://doi.org/10.1038/nature12650' target='_blank'>View</a>&nbsp;&nbsp;&nbsp;
								<a href='https://www.ncbi.nlm.nih.gov/pubmed/24185008' target='_blank'>Pubmed</a>
							</li>
						</ul>
					")
				)
			),
			
	##~~~~~~~~~~~~~~~~~~~~~~~~~~  COVERAGE  ~~~~~~~~~~~~~~~~~~~~~~~~~~##
			
			
			conditionalPanel(condition="input.tabs=='Coverage'",
				div(id='coverageSide',
					div(id='loadBedDiv',
						class="inline",
						style="padding-bottom:30px;width:75%;margin-right:10px;",
						HTML("<div class='myTooltip'><label>Load a bed file:</label><span class='myTooltiptext'>
							File must be produced with the CANmapper script.
							</span></div>"),
						div(class="inline",
							fileInput("bedFile", NULL,multiple=F,accept=".bed",buttonLabel = "Browse...", placeholder="No file selected")
						)
					),
					div(id='analyseOrReset',
						class="inline",
						style="padding-top:25px;"
					),
					div(id='exampleCoverageDiv',
						style="margin-top:-10px;",
						HTML("<div style='padding-bottom:5px;'><b>Or</b></div>"),
						actionButton('exampleCoverage',"Load example")
					)
				)
			),

	##~~~~~~~~~~~~~~~~~~~~~~~~~~~  RATIO  ~~~~~~~~~~~~~~~~~~~~~~~~~~~##

			conditionalPanel(condition="input.tabs=='Ratio'",
				div(id='ratioSide',
					div(id='loadRatioDiv',
						style="margin-bottom:30px;width:75%;margin-right:10px;",
						HTML("<div class='myTooltip'><label>Load a saved ratio file:</label><span class='myTooltiptext'>
							File must be produced using this page
							</span></div>"),
						div(class='inline',
							fileInput("ratioFile",NULL,multiple=F,buttonLabel = "Browse...", placeholder="No file selected")
						)
					),
					div(id='exampleRatioDiv',
						style="margin-top:-10px;",
						HTML("<div style='padding-bottom:5px;'><b>Or</b></div>"),
						actionButton('exampleRatio',"Load example")
					)
				)
			),

	##~~~~~~~~~~~~~~~~~~~~~~~~~~~  PLOT  ~~~~~~~~~~~~~~~~~~~~~~~~~~~##

			conditionalPanel(condition="input.tabs=='Plot'",
				div(id='plotSide',
					div(
						id='loadPlotDiv',
						style="margin-bottom:30px;width:75%;margin-right:10px;",
						HTML("<div class='myTooltip'><label>Load a saved ratios file:</label><span class='myTooltiptext'>
							File must be produced using this page</span></div>"
						),
						fileInput('plotFile',NULL,multiple=F,buttonLabel = "Browse...", placeholder="No file selected")
					),
					div(
						id='examplePlotDiv',
						style="margin-top:-10px;",
						HTML("<div style='padding-bottom:5px;'><b>Or</b></div>"),
						actionButton('examplePlot',"Load example")
					),
					div(
						id='plotSideCtrls',
						div(id='samples',style='padding-bottom:15px;')
					)
				)
			),
			
			
			
	##~~~~~~~~~~~~~~~~~~~~~~~~~~~  STATS  ~~~~~~~~~~~~~~~~~~~~~~~~~~~##

			conditionalPanel(condition="input.tabs=='Stats'",
				div(id="statsSide",
					div(
						id='exampleStatsDiv',
						style='padding-top:30px;',
						actionButton('exampleStats',"Load example")
					),
					div(id='statSamples',style='padding-bottom:15px;')
				)
			),
			width=3
		),


		mainPanel(
			tabsetPanel(
				tabPanel("About",
					div(id="aboutMain",
						style="padding:20px;",
						HTML('
						<div class="col-sm-5"><h3>Replication profiling</h3>
						<p style="padding-top:10px;text-align:justify;color:#404040;">
							At the population level, DNA replication follows a defined programme that is subject to changes in differentiation 
							and disease. However, at the single cell level, the process is stochastic. Early and efficient origins of 
							replication fire at the onset of the synthetic (S) phase, initiating the process of DNA replication. Origin-less 
							regions are passively replicated late in S phase by forks originating from neighbouring origins. By combining FACS 
							with high-throughput sequencing, it is possible to capture the average replication dynamics of a cell 
							population. The protocol includes fixation of asynchronous cells, FACS enrichment of the S phase population, as 
							well as non-replicating control cells, deep sequencing of the extracted DNA, and computational analysis. Therefore, 
							it is a direct, copy number-based quantitative approach for investigating the replication dynamics of the whole 
							nuclear genome without perturbing the cells. It is applicable to any single-celled model organism with a stable 
							genome that can be sorted by FACS. Currently, it is mostly used on organisms with small genomes. However, 
							up-and-coming advancement of deep sequencing approaches will facilitate the analysis of organisms with genomes 
							of any size. This is of broad interest not only to DNA replication, but also the fields of DNA repair, evolution 
							research, pre-clinical studies of anti-cancer drugs and the characterisation of novel organisms. Our comprehensive 
							pipeline includes quality controls before the deep sequencing, normalisation between samples using a foreign DNA 
							spike-in approach, and computational scripts to analyse and plot deep sequencing data, as well as a statistical 
							test to assess differences between replication profiles. 
						</p></div>
						<div class="col-sm-7">
							<div style="overflow:hidden;max-width:850px;min-width:500px;padding:50px 10px;"><img src="outline.jpg" width="100%"></div>
						</div>')
					)
				),
				tabPanel("Coverage",
					div(id="coverageMain",
						HTML("
							<div id='coverageDescription' style='padding:50px;text-align:justify;color:#404040;'>
								<p>Use the menu on the left to either load example data or upload your own. The uploaded bed file may or may not have a 
								header and <b>must</b> contain 5 columns (<b><i>chrom, chromStart, chromEnd, name, score</i></b>) where the first three columns 
								define the genomic bin and the \"score\" column is used for storing coverage data (reads per bin). The \"name\" 
								column will be repurposed to store the bed file's name, therefore bed file names should be descriptive.</p>
								<p>When the file is uploaded, a snippet of its modified content is displayed initially - please make sure everything looks right
								before continuing.</p>
								<p><b>Tip:</b> Hover over headings in the control panel to discover tooltips!</p>
							</div>
						")
					)
				),
				tabPanel("Ratio",
					div(id="ratioMain",
						HTML("
							<div id='ratioDescription' style='padding:50px;text-align:justify;color:#404040;'>
								<p>Use the menu on the left to either load example data or upload your own (created earlier using this page!). 
								Once you have saved at least one sample for each replicating and non-replicating sample type, this page 
								will display elements for making a new ratio. Initially, the ratio is normalised by the total read 
								number and will have a distribution around one.</p>
								<p>In the case of full range S phase samples (sorted whole S phase or synchronised 
								S phase population, where at least some regions are completely replicated), <b>automatic normalisation</b> may be used. 
								This scales the data to lie between one and two, based on minimising the sum of data points outside of this region. 
								<b>Trimming</b> should be done if there are ratio values far outside of the main population, as they will skew the automatic 
								normalisation. A range of 0.5-1.5 is a very safe starting point.</p>
								<p>If replicating samples come from early S phase timepoints of a cell cycle experiment, or an asynchronous cell 
								culture (marker frequency analysis), <b>manual normalisation</b> should be used. For example, if the asynchronous population contains 20% 
								of cells in S phase, a factor of 1.2 should be used. Likewise, if cells from a synchronised S phase population are 
								15% replicated, a factor of 1.15 should be used.</p>
							</div>
						")
					)
				),
				tabPanel("Plot",
					div(id="plotMain",
						HTML("
							<div id='plotDescription' style='padding:50px;text-align:justify;color:#404040;'>
							
								<p>Use menu on the left to either load an example data or upload your own multiple ratios file (created earlier 
								using this page!). Any ratios saved in the previous tab should appear on the left, along with some controls.</p>
								
								<p><b>Additional features</b> may be added to the plot. Typically, replication origins are plotted as circles and centromere 
								location as vertical lines. Rectangles can be used to highlight a big region of a chromosome, while pointers can be used to 
								pinpoint a specific locus. Use the \"name\" field in the bed file to name the feature.</p>
								
								<p> Run <b>smoothing</b> to plot smoothed data. Due to the discrete nature of the data (separate chromosomes, 
								missing bins), smoothing is done in groups. <label>Group size</label> controls a minimum number of datapoints 
								to smooth (each group must contain a minimum of 4 datapoints for the spline algorithm to work). 
								<label>Split</label> value controls the number of missing bins along a chromosome to initiate a new 
								group.</p>
								
								<p>Use <b>plotting controls</b> to zoom into a particular chromosome (or a subchromosomal region) or
								choose different plot type (note that the only suitable plot type to display both raw and smooth data is scatter plot). 
								Y axis limits allow to plot marker frequency analysis data.</p>
							</div>
						")
					)
				),
				tabPanel("Stats",
					div(id="statsMain",
						HTML("
							<div id='statsDescription' style='padding:50px;text-align:justify;color:#404040;'>
								<p>Use menu on the left to load example data.</p>
								<p>Any ratios saved in the 'Ratio' tab should appear on the left (when at least two are available). Plotting controls are 
								similar to the 'Plot' tab.</p>
								<p>The statistical analysis is based on <a href='https://en.wikipedia.org/wiki/Standard_score' target='_blank'>z-scores</a>. 
								Statistically significant differences between two replication profiles are separated into bins of 0.99-0.999 significance 
								(p-value**) and above 0.999 significance (p-value***).</p>
							</div>
						")
					)
				),id="tabs"
			), width=9
		)
	),
	HTML("<div class='footer'><span style='height:50px;line-height:50px;vertical-align:middle;'>
		<a href='http://nieduszynski.org/' target='_blank'>Nieduszynski lab</a></span></div>")
)
