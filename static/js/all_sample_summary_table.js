d3.csv("report_data/sample_summary.csv", function (samples) {

    // original idea stolen from http://bl.ocks.org/gka/17ee676dc59aa752b4e6
    
    // column definitions
    var columns = [
        { head: 'Sample name', cl: 'title', csv_title: 'name' },
        //{ head: 'thing2', cl: 'title', csv_title: 'thing' }
    ];

    var table = d3.select('body')
        .append('table')
    
    // create table header
    table.append('thead').append('tr')
        .selectAll('th')
        .data(columns).enter()
        .append('th')
        .attr('class', function (c) { return c.cl })
        .text( function (c) { return c.head } );

    // create table body
    table.append('tbody')
        .selectAll('tr')
        .data(samples).enter()
        .append('tr')
        .selectAll('td')
        .data(function(row, i) {
            return columns.map(function(c) {
                // compute cell values for this specific row
                
                var cell = {};
                d3.keys(c).forEach( function(k) {
                    cell['html'] = row[c.csv_title]
                });
                return cell;
            });
        }).enter()
        .append('td')
        .html(function (s) { return s.html })
})