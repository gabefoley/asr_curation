

dataset_ranges_text = f"""
    DATASET_RANGE
#Colored/labeled range datasets allow the highlighting of various clades or leaf ranges by using colored boxes or brackets.

#lines starting with a hash sign are comments and ignored during parsing
#=================================================================#
#                    MANDATORY SETTINGS                           #
#=================================================================#
#select the separator which is used to delimit the data below (TAB,SPACE or COMMA).This separator must be used throughout this file.
#SEPARATOR TAB
#SEPARATOR SPACE
SEPARATOR COMMA

#label is used in the legend table (can be changed later)
DATASET_LABEL, <custom_dataset_label>_ranges

#dataset color in the legend table
COLOR,#ffff00


#=================================================================#
#                    OPTIONAL SETTINGS                            #
#=================================================================#

#=================================================================#
#        all optional settings can be set or changed later        #
#           in the web interface (under 'Datasets' tab)           #
#=================================================================#

#RANGE_TYPE defines how the rages will be visualized:
   #box: standard colored box/polygon. Various LINE_? fields in the range definition will be used for the border style. 
   #bracket: a line or bracket outside the tree

RANGE_TYPE,box
#
#Box/polygon specific options, used when RANGE_TYPE is 'box'
#

#specify what the range boxes will cover: 'label','clade' or 'tree'
RANGE_COVER,clade

#simplify or smooth polygons when in unrooted display mode: 'none', 'simplify' or 'smooth'
UNROOTED_SMOOTH,simplify

#when RANGE_COVER is set to 'clade' or 'tree', you can disable the covering of labels (ie. limiting the boxes to the tree structure only)
COVER_LABELS,0

#if set to 1, ranges will cover any displayed extrernal datasets as well
COVER_DATASETS,0

#if set to 1, size of the boxes will be extended to fit their labels
FIT_LABELS,0

#
#Bracket specific options, used when RANGE_TYPE is 'bracket'
#

#bracket style can be: 'none','square' or 'curved'
BRACKET_STYLE,square

#size of the bracket ends (for 'square' or 'curved' brackets)
BRACKET_SIZE,20

#shift the bracket position horizontally
BRACKET_SHIFT,50

#if set to 1, brackets will be displayed behind the last visible external dataset
BRACKET_BEHIND_DATASETS,1

#
#Options related to range labels
#

SHOW_LABELS,1

#the position of the label in the range box (or relative to the bracket): 'top-left','top-center','top-right',
#                                                                         'center-left','center-center','center-right',
#                                                                         'bottom-left','bottom-center','bottom-right'
LABEL_POSITION,bottom-right

#Display the labels vertically. In circular display mode (or with brackets in unrooted display mode), labels will be aligned to the circle
LABELS_VERTICAL,0

#labels remain straight, regardless of the tree rotation or other rotation parameters
STRAIGHT_LABELS,0

#rotate all labels by the specified angle
LABEL_ROTATION,0

#shift all labels horizontally and/or vertically
LABEL_SHIFT_X,0
LABEL_SHIFT_Y,0

#add a colored outline to the label font; useful when displaying labels over similarly colored boxes (e.g. black font on a dark box)
LABEL_OUTLINE_WIDTH,0
LABEL_OUTLINE_COLOR,#ffffff

#multiply the size of all labels by this factor
LABEL_SIZE_FACTOR,1

#shrink range boxes or brackets vertically, to introduce spacing between neighbouring ranges
VERTICAL_SHRINK,0

#Each dataset can have a legend, which is defined using LEGEND_XXX fields below
#For each row in the legend, there should be one shape, color and label.
#Optionally, you can define an exact legend position using LEGEND_POSITION_X and LEGEND_POSITION_Y. To use automatic legend positioning, do NOT define these values
#Optionally, shape scaling can be present (LEGEND_SHAPE_SCALES). For each shape, you can define a scaling factor between 0 and 1.
#To order legend entries horizontally instead of vertically, set LEGEND_HORIZONTAL to 1
#Shape should be a number between 1 and 6, or any protein domain shape definition.
#1: square
#2: circle
#3: star
#4: right pointing triangle
#5: left pointing triangle
#6: checkmark

#LEGEND_TITLE,Dataset legend
#LEGEND_POSITION_X,100
#LEGEND_POSITION_Y,100
#LEGEND_HORIZONTAL,0
#LEGEND_SHAPES,1,2,3
#LEGEND_COLORS,#ff0000,#00ff00,#0000ff
#LEGEND_LABELS,value1,value2,value3
#LEGEND_SHAPE_SCALES,1,1,0.5



#Internal tree nodes can be specified by using IDs directly, or by using the 'last common ancestor' method described in iTOL help pages
#=================================================================#
#       Actual data follows after the "DATA" keyword              #
#=================================================================#
#the following fields are available in each line:

#START_NODE_ID,END_NODE_ID,FILL_COLOR,GRADIENT_COLOR,LINE_COLOR,LINE_STYLE,LINE_WIDTH,LABEL_TEXT,LABEL_COLOR,LABEL_SIZE_FACTOR,LABEL_STYLE

#The range is defined through START_NODE_ID and END_NODE_ID.
#If GRADIENT_FILL color is defined, the box will be filled with a gradient from FILL_COLOR to GRADIENT_COLOR.  Brackets will also be visualized as gradients.
#LINE_COLOR will be used for the box/polygon border, or for the brackets. If not specified, FILL_COLOR will be used instead
#LINE_STYLE can be 'solid', 'dashed' or 'dotted'
#LABEL_STYLE can be 'normal', 'bold', 'italic' or 'bold-italic'

DATA
#Examples
#a range between leaves 9606 and 184922, filled with a gradient from white (#ffffff) to red (#ff0000), with a 2px dashed black (#000000) border and a blue (#0000ff) italic label
#9606,184922,#ffffff,#ff0000,#000000,dashed,2,Example range,#0000ff,1,italic
    """

dataset_colorstrip_text = f"""DATASET_COLORSTRIP
#In colored strip datasets, each ID is associated to a color box/strip and can have an optional label. Color can be specified in hexadecimal, RGB or RGBA notation. When using RGB or RGBA notation, you cannot use COMMA as the dataset separator

#lines starting with a hash are comments and ignored during parsing

#=================================================================#
#                    MANDATORY SETTINGS                           #
#=================================================================#
#select the separator which is used to delimit the data below (TAB,SPACE or COMMA).This separator must be used throughout this file.

#SEPARATOR TAB
SEPARATOR COMMA
#SEPARATOR SPACE

#label is used in the legend table (can be changed later)
DATASET_LABEL,<custom_dataset_label>_colorstrip

#dataset color (can be changed later)
COLOR,#ff0000

#=================================================================#
#                    OPTIONAL SETTINGS                            #
#=================================================================#

#If COLOR_BRANCHES is set to 1, branches of the tree will be colored according to the colors of the strips above the leaves.
#When all children of a node have the same color, it will be colored the same, ie. the color will propagate inwards towards the root.
COLOR_BRANCHES,0


#=================================================================#
#     all other optional settings can be set or changed later     #
#           in the web interface (under 'Datasets' tab)           #
#=================================================================#

#Each dataset can have a legend, which is defined using LEGEND_XXX fields below
#For each row in the legend, there should be one shape, color and label.
#Optionally, you can define an exact legend position using LEGEND_POSITION_X and LEGEND_POSITION_Y. To use automatic legend positioning, do NOT define these values
#Optionally, shape scaling can be present (LEGEND_SHAPE_SCALES). For each shape, you can define a scaling factor between 0 and 1.
#To order legend entries horizontally instead of vertically, set LEGEND_HORIZONTAL to 1
#Shape should be a number between 1 and 6, or any protein domain shape definition.
#1: square
#2: circle
#3: star
#4: right pointing triangle
#5: left pointing triangle
#6: checkmark

#LEGEND_TITLE Dataset_legend
#LEGEND_POSITION_X 100
#LEGEND_POSITION_Y 100
#LEGEND_HORIZONTAL 0
#LEGEND_SHAPES 1 1 2 2
#LEGEND_COLORS #ff0000 #00ff00 rgba(0,255,0,0.5) #0000ff
#LEGEND_LABELS value1 value2 value3 value4
#LEGEND_SHAPE_SCALES 1 1 0.5 1

#width of the colored strip
#STRIP_WIDTH 25

#left margin, used to increase/decrease the spacing to the next dataset. Can be negative, causing datasets to overlap.
#MARGIN 0

#border width; if set above 0, a border of specified width (in pixels) will be drawn around the color strip 
#BORDER_WIDTH 0

#border color; used when BORDER_WIDTH is above 0
#BORDER_COLOR #0000ff

#if set to 1, border will be drawn completely around each colored strip box
#COMPLETE_BORDER 0

#always show internal values; if set, values associated to internal nodes will be displayed even if these nodes are not collapsed. It could cause overlapping in the dataset display.
#SHOW_INTERNAL 0


#display or hide the individual label inside each colored strip (when defined in the data below)
#SHOW_STRIP_LABELS 1

#position of the strip label within the box; 'top', 'center' or 'bottom'
#STRIP_LABEL_POSITION center

#strip label size factor (relative to the tree leaf labels)
#STRIP_LABEL_SIZE_FACTOR 1


#rotation of the strip labels; used only in rectangular tree display mode
#STRIP_LABEL_ROTATION 0

#strip label shift in pixels (positive or negative)
#STRIP_LABEL_SHIFT 0

#STRIP_LABEL_COLOR #000000

#draw a black outline around the text (width in pixels)
#STRIP_LABEL_OUTLINE 0.5

#display or hide the dataset label above the colored strip
#SHOW_LABELS 1

#dataset label size factor
#SIZE_FACTOR 1

#dataset label rotation
#LABEL_ROTATION 0

#dataset label shift in pixels (positive or negative)
#LABEL_SHIFT 0

#Internal tree nodes can be specified using IDs directly, or using the 'last common ancestor' method described in iTOL help pages

#=================================================================#
#       Actual data follows after the "DATA" keyword              #
#=================================================================#
DATA

#Examples:
#assign a red colored strip to leaf 9606, with label 'Human'
#9606 #ff0000 Human

#assign a green, semi-transparent (alpha 0.5) strip to an internal node, without any label. If 'Show internal values' is set to 'No', this will only be displayed if the node is collapsed. 
#9606|5664 rgba(0,255,0,0.5)
"""

dataset_text_text = f"""
DATASET_TEXT
#In text datasets, each ID is associated to text label, which can be displayed directly on the node branch, or outside the tree
#lines starting with a hash are comments and ignored during parsing
#=================================================================#
#                    MANDATORY SETTINGS                           #
#=================================================================#
#select the separator which is used to delimit the data below (TAB,SPACE or COMMA).This separator must be used throughout this file.
#SEPARATOR TAB
#SEPARATOR SPACE
SEPARATOR COMMA

#label is used in the legend table (can be changed later)
DATASET_LABEL,<custom_dataset_label>_text

#dataset color (can be changed later)
COLOR,#ff0000

#=================================================================#
#                    OPTIONAL SETTINGS                            #
#=================================================================#

#=================================================================#
#     all other optional settings can be set or changed later     #
#           in the web interface (under 'Datasets' tab)           #
#=================================================================#

#left margin, used to increase/decrease the spacing to the next dataset. Can be negative, causing datasets to overlap. Used only for text labels which are displayed on the outside
MARGIN,0

#applies to external text labels only; if set, text labels associated to internal nodes will be displayed even if these nodes are not collapsed. It could cause overlapping in the dataset display.
SHOW_INTERNAL,0

#Rotate all labels by the specified angle
ALL_LABELS_ROTATION,0

#By default, internal labels will be placed above the branches. If LABELS_BELOW is set to 1, labels will be below the branches
LABELS_BELOW,1

#Shift internal labels vertically by this amount of pixels (positive or negative)
VERTICAL_SHIFT,-20

#If set to 1, tree rotation will not influence the individual label rotation
STRAIGHT_LABELS,0

#applies to external text labels only; If set to 1, labels will be displayed in arcs aligned to the tree (in circular mode) or vertically (in normal mode). All rotation parameters (global or individual) will be ignored.
ALIGN_TO_TREE,0

#font size factor; For external text labels, default font size will be slightly less than the available space between leaves, but you can set a multiplication factor here to increase/decrease it (values from 0 to 1 will decrease it, values above 1 will increase it)
SIZE_FACTOR,1

#add extra horizontal shift to the external labels. Useful in unrooted display mode to shift text labels further away from the node labels.
EXTERNAL_LABEL_SHIFT,0


#display or hide the dataset label above the external labels column
#SHOW_LABELS 1

#dataset label size factor
#LABEL_SIZE_FACTOR 1

#dataset label rotation
#LABEL_ROTATION 0

#dataset label shift in pixels (positive or negative)
#LABEL_SHIFT 0


#Internal tree nodes can be specified using IDs directly, or using the 'last common ancestor' method described in iTOL help pages
#=================================================================#
#       Actual data follows after the "DATA" keyword              #
#=================================================================#
#the following fields are possible for each node:
#ID,label,position,color,style,size_factor,rotation

#position defines the position of the text label on the tree:
#  -1 = external label
#  a number between 0 and 1 = internal label positioned at the specified value along the node branch (for example, position 0 is exactly at the start of node branch, position 0.5 is in the middle, and position 1 is at the end)
#style can be 'normal',''bold','italic' or 'bold-italic'
#size factor will be multiplied with the standard font size

DATA
#Examples

#node 9598 will have an external label 'Pan troglodytes' in bold red and twice the size of standard labels
#9598,Pan troglodytes,-1,#ff0000,bold,2,0

#node 9606 will have an external label with multiple mixed styles
#9606,<bi color='#006600'>Homo </bi><i>sapiens</i><sup size='0.5' color='#999999'>citation</sup>,-1,#000000,normal,1,0

#node 4530 will have an internal label 'Oryza sativa' in bold italic blue, starting directly over the node
#4530,Oryza sativa,0,#0000ff,bold-italic,1,0
"""



lab_assays_text = f"""
    DATASET_EXTERNALSHAPE
#Nodes have multiple values associated with them. Values will be displayed as geometric shapes of different sizes in columns outside the tree.
#Highest value in the dataset will have the largest size, and all others will be scaled down proportionally.
#lines starting with a hash are comments and ignored during parsing

#=================================================================#
#                    MANDATORY SETTINGS                           #
#=================================================================#
#select the separator which is used to delimit the data below (TAB,SPACE or COMMA).This separator must be used throughout this file.
#SEPARATOR TAB
#SEPARATOR SPACE
SEPARATOR COMMA

#label is used in the legend table (can be changed later)
DATASET_LABEL, Lab Assays

#dataset color (can be changed later)
COLOR,#ff0000

#define colors for each individual field column (use hexadecimal, RGB or RGBA notation; if using RGB/RGBA, COMMA cannot be used as SEPARATOR)
FIELD_COLORS,#ff0000,#ff0000,#ff0000,#00ff00,#00ff00,#00ff00,#0000ff,#0000ff,#0000ff,#800080,#800080,#800080,#ff6666,#ff6666,#ff6666,#66ff99,#66ff99,#66ff99,#6666ff,#6666ff,#6666ff,#ffcc99,#ffcc99,#ffcc99,#ff9966,#ff9966,#ff9966,#99ff99,#99ff99,#99ff99,#9999ff,#9999ff,#9999ff,#ff9966,#ff9966,#ff9966

#field labels
FIELD_LABELS,Activity_Zn_PnP, Activity_Mg_PnP, Activity_Mn_PnP, Activity_Zn_BnP, Activity_Mg_BnP, Activity_Mn_BnP, Activity_Zn_TnP, Activity_Mg_TnP, Activity_Mn_TnP, Activity_Zn_4NPC, Activity_Mg_4NPC, Activity_Mn_4NPC, Activity_Zn_4NPS, Activity_Mg_4NPS, Activity_Mn_4NPS, Activity_Zn_4NPA, Activity_Mg_4NPA, Activity_Mn_4NPA, Activity_Zn_4NPB, Activity_Mg_4NPB, Activity_Mn_4NPB, Activity_Zn_4NPH, Activity_Mg_4NPH, Activity_Mn_4NPH, Activity_Zn_DVL, Activity_Mg_DVL, Activity_Mn_DVL, Activity_Zn_DHM, Activity_Mg_DHM, Activity_Mn_DHM, Activity_Zn_Nitrocefin, Activity_Mg_Nitrocefin, Activity_Mn_Nitrocefin, Activity_Zn_SLG, Activity_Mg_SLG, Activity_Mn_SLG

#=================================================================#
#                    OPTIONAL SETTINGS                            #
#=================================================================#


#=================================================================#
#     all other optional settings can be set or changed later     #
#           in the web interface (under 'Datasets' tab)           #
#=================================================================#

#Each dataset can have a legend, which is defined using LEGEND_XXX fields below
#For each row in the legend, there should be one shape, color and label.
#Optionally, you can define an exact legend position using LEGEND_POSITION_X and LEGEND_POSITION_Y. To use automatic legend positioning, do NOT define these values
#Optionally, shape scaling can be present (LEGEND_SHAPE_SCALES). For each shape, you can define a scaling factor between 0 and 1.
#To order legend entries horizontally instead of vertically, set LEGEND_HORIZONTAL to 1
#Shape should be a number between 1 and 6, or any protein domain shape definition.
#1: square
#2: circle
#3: star
#4: right pointing triangle
#5: left pointing triangle
#6: checkmark

LEGEND_TITLE,Total lab assays
LEGEND_SHAPES,3,3,3
LEGEND_COLORS,#ff0000,#0000ff,#00ff00
LEGEND_LABELS,Completed - <completed>, Total - <total>, Percentage completed - <percentage>
#LEGEND_SHAPE_SCALES,1,0.5,1,1,1
#LEGEND_SHAPE_INVERT,1,0,0,0,0
#left margin, used to increase/decrease the spacing to the next dataset. Can be negative, causing datasets to overlap.
#MARGIN,0

#left margin, used to increase/decrease the spacing to the next dataset. Can be negative, causing datasets to overlap.
#MARGIN,0

#always show internal values; if set, values associated to internal nodes will be displayed even if these nodes are not collapsed. It could cause overlapping in the dataset display.
#SHOW_INTERNAL,0

#show dashed lines between leaf labels and the dataset
DASHED_LINES,1

#shape height factor; Default shape height will be slightly less than the available space between leaves, but you can set a multiplication factor here to increase/decrease it (values from 0 to 1 will decrease it, values above 1 will increase it)
#HEIGHT_FACTOR,1

#vertical and horizontal grids can be displayed connecting the shapes
HORIZONTAL_GRID,1
VERTICAL_GRID,1

#horizontal spacing between shape columns
SHAPE_SPACING,5

#Shape types:
#1: square
#2: circle
#3: star
#4: right pointing triangle
#5: left pointing triangle
SHAPE_TYPE,2,2,2,1,1,1,2,2,2

#if set to 0, only outlines will be shown
#COLOR_FILL,1

#if set to 1, actual numeric value will be show in the center of each shape
SHOW_VALUES,1

#display or hide the text labels above each field column
SHOW_LABELS,1

#text label size factor
#SIZE_FACTOR,1

#text label rotation
#LABEL_ROTATION,0

#text label shift in pixels (positive or negative)
#LABEL_SHIFT,0


#=================================================================#
#       Actual data follows after the "DATA" keyword              #
#=================================================================#
DATA

"""

dataset_style_text = f"""
DATASET_STYLE
#Style datasets allow the customization of branch and leaf label colors and styles.

#lines starting with a hash are comments and ignored during parsing
#=================================================================#
#                    MANDATORY SETTINGS                           #
#=================================================================#
#select the separator which is used to delimit the data below (TAB,SPACE or COMMA).This separator must be used throughout this file.
#SEPARATOR TAB
#SEPARATOR SPACE
SEPARATOR COMMA

#label is used in the legend table (can be changed later)
DATASET_LABEL, <custom_dataset_label>_dataset_style

#dataset color (can be changed later)
COLOR,#ffff00

#=================================================================#
#                    OPTIONAL SETTINGS                            #
#=================================================================#

#Each dataset can have a legend, which is defined using LEGEND_XXX fields below
#For each row in the legend, there should be one shape, color and label.
#Optionally, you can define an exact legend position using LEGEND_POSITION_X and LEGEND_POSITION_Y. To use automatic legend positioning, do NOT define these values
#Optionally, shape scaling can be present (LEGEND_SHAPE_SCALES). For each shape, you can define a scaling factor between 0 and 1.
#To order legend entries horizontally instead of vertically, set LEGEND_HORIZONTAL to 1
#Shape should be a number between 1 and 6, or any protein domain shape definition.
#1: square
#2: circle
#3: star
#4: right pointing triangle
#5: left pointing triangle
#6: checkmark

#LEGEND_TITLE,Dataset legend
#LEGEND_POSITION_X,100
#LEGEND_POSITION_Y,100
#LEGEND_HORIZONTAL,0
#LEGEND_SHAPES,1,2,3
#LEGEND_COLORS,#ff0000,#00ff00,#0000ff
#LEGEND_LABELS,value1,value2,value3
#LEGEND_SHAPE_SCALES,1,1,0.5

#Internal tree nodes can be specified using IDs directly, or using the 'last common ancestor' method described in iTOL help pages
#=================================================================#
#       Actual data follows after the "DATA" keyword              #
#=================================================================#
#the following fields are required for each node:
#ID,TYPE,WHAT,COLOR,WIDTH_OR_SIZE_FACTOR,STYLE,BACKGROUND_COLOR

# TYPE: can be either 'branch' or 'label'. 'branch' will apply customizations to the tree branches, while 'labels' apply to the leaf text labels
# WHAT: can be either 'node' or 'clade', only relevant for internal tree nodes. 'Node' will apply the customization only to a single node, while 'clade' will apply to all child nodes as well.
# COLOR: can be in hexadecimal, RGB or RGBA notation. If RGB or RGBA are used, dataset SEPARATOR cannot be comma.
# WIDTH_OR_SIZE_FACTOR: for type 'branch', specifies the relative branch width, compared to the global branch width setting.
#                       for type 'label', specifies the relative font size, compared to the global font size
# STYLE: for type 'branch', can be either 'normal' or 'dashed'
#        for type 'label', can be one of 'normal', 'bold', 'italic' or 'bold-italic'
# BACKGROUND_COLOR (optional): only relevant for type 'label', specifies the color of the label background. The value is optional.


DATA

#Examples

#a single internal node's branch will be colored red with double branch width and dashed line
#9606|184922,branch,node,#ff0000,2,dashed

#node 9606 will have its label displayed in blue with bold italic font, and with yellow background
#9606,label,node,#0000ff,1,bold-italic,#ffff00

#a clade starting at internal node 2190|2287 will have all its branches colored green
#2190|2287,branch,clade,#00ff00,1,normal

#all leaf labels in a clade will be displayed in red
#2097|1502,label,clade,#ff0000,1,normal
"""

shape_text = '''
DATASET_EXTERNALSHAPE
#Nodes have multiple values associated with them. Values will be displayed as geometric shapes of different sizes in columns outside the tree.
#Highest value in the dataset will have the largest size, and all others will be scaled down proportionally.
#lines starting with a hash are comments and ignored during parsing

#=================================================================#
#                    MANDATORY SETTINGS                           #
#=================================================================#
#select the separator which is used to delimit the data below (TAB,SPACE or COMMA).This separator must be used throughout this file.
#SEPARATOR TAB
#SEPARATOR SPACE
SEPARATOR COMMA

#label is used in the legend table (can be changed later)
DATASET_LABEL, <custom_dataset_label>_shapes

#dataset color (can be changed later)
COLOR,red

#define colors for each individual field column (use hexadecimal, RGB or RGBA notation; if using RGB/RGBA, COMMA cannot be used as SEPARATOR)
FIELD_COLORS,<custom_color>,#00ff00,#0000ff

#field labels
FIELD_LABELS,<col>

#=================================================================#
#                    OPTIONAL SETTINGS                            #
#=================================================================#


#=================================================================#
#     all other optional settings can be set or changed later     #
#           in the web interface (under 'Datasets' tab)           #
#=================================================================#

#Each dataset can have a legend, which is defined using LEGEND_XXX fields below
#For each row in the legend, there should be one shape, color and label.
#Optionally, you can define an exact legend position using LEGEND_POSITION_X and LEGEND_POSITION_Y. To use automatic legend positioning, do NOT define these values
#Optionally, shape scaling can be present (LEGEND_SHAPE_SCALES). For each shape, you can define a scaling factor between 0 and 1.
#To order legend entries horizontally instead of vertically, set LEGEND_HORIZONTAL to 1
#Shape should be a number between 1 and 6, or any protein domain shape definition.
#1: square
#2: circle
#3: star
#4: right pointing triangle
#5: left pointing triangle
#6: checkmark

#LEGEND_TITLE,Dataset legend
#LEGEND_POSITION_X,100
#LEGEND_POSITION_Y,100
#LEGEND_HORIZONTAL,0
#LEGEND_SHAPES,1,2,3
#LEGEND_COLORS,#ff0000,#00ff00,#0000ff
#LEGEND_LABELS,value1,value2,value3
#LEGEND_SHAPE_SCALES,1,1,0.5

#left margin, used to increase/decrease the spacing to the next dataset. Can be negative, causing datasets to overlap.
#MARGIN,0

#always show internal values; if set, values associated to internal nodes will be displayed even if these nodes are not collapsed. It could cause overlapping in the dataset display.
#SHOW_INTERNAL,0

#show dashed lines between leaf labels and the dataset
DASHED_LINES,1

#shape height factor; Default shape height will be slightly less than the available space between leaves, but you can set a multiplication factor here to increase/decrease it (values from 0 to 1 will decrease it, values above 1 will increase it)
#HEIGHT_FACTOR,1

#vertical and horizontal grids can be displayed connecting the shapes
#HORIZONTAL_GRID,1
#VERTICAL_GRID,1

#horizontal spacing between shape columns
#SHAPE_SPACING,10

#Shape types:
#1: square
#2: circle
#3: star
#4: right pointing triangle
#5: left pointing triangle
#SHAPE_TYPE,2

#if set to 0, only outlines will be shown
#COLOR_FILL,1

#if set to 1, actual numeric value will be show in the center of each shape
#SHOW_VALUES,1

#display or hide the text labels above each field column
SHOW_LABELS,1

#text label size factor
#SIZE_FACTOR,1

#text label rotation
#LABEL_ROTATION,0

#text label shift in pixels (positive or negative)
#LABEL_SHIFT,0


#=================================================================#
#       Actual data follows after the "DATA" keyword              #
#=================================================================#
DATA
#Internal tree nodes can be specified using IDs directly, or using the 'last common ancestor' method described in iTOL help pages
#ID1,value1,value2,value3
#ID2,value4,value5,value6
#9606,10,10,20,40
#LEAF1|LEAF2,50,60,80,90
'''