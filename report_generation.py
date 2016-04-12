def generate_index(analysispipe):
    
    outputStr = """
    <!DOCTYPE HTML>
        <html lang="en-US">
            <head>
                <meta charset="UTF-8">
                <script src="http://d3js.org/d3.v3.js"></script>
                <link rel="stylesheet" href="static/css/style.css">
            </head>
            <body>
            <script src="https://rawgit.com/gka/d3-jetpack/master/d3-jetpack.js"></script>
            <script src="/static/js/all_sample_summary_table.js"></script>
            <div id="all_sample_summary_table"></div>"""
    outputStr += """
            </body>
        </html>"""
        
    return outputStr