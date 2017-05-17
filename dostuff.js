$.urlParam = function(name){
              var results = new RegExp('[\?&]' + name + '=([^&#]*)').exec(window.location.href);
              return results[1] || 0;
          }

          $(document).ready(function () {

            if ($.urlParam('axisTypeOptions') == "minorAxis") {

              $("#minorAxisForm").show();

            } else if($.urlParam('axisTypeOptions') == "majorAxis") {

              $("#majorAxisForm").show();
            }

          $("#majorAxisForm").submit(function( event ) {
            var $inputs = $('#majorAxisForm :input');
            var values = [];
            $inputs.each(function() {
              values.push(parseFloat($(this).val()));
            });

            console.log(values);
			var result = maincode(values[0], values[1], values[2], values[3], values[4], values[5],0, values[6], values[7], values[8]);
            console.log(result);
			console.log(result[0].length);
			plot_chart(result);
			/*var json_result = JSON.stringify(result);
			createCookie('mycookie', json_result);
			window.location.href = "chart.html";*/
			event.preventDefault();
          });

          $("#minorAxisForm").submit(function( event ) {
            var $inputs = $('#minorAxisForm :input');
            var values = [];
            $inputs.each(function() {
              values.push(parseFloat($(this).val()));
            });

            console.log(values);
			var result = maincodeminor(values[0], values[1], values[2], values[3], values[4], values[5], 1, values[6], values[7], values[8]);
            console.log(result);
			console.log(result[0].length);
			plot_chart(result);
			/*var json_result = JSON.stringify(result);
			createCookie('mycookie', json_result);
			window.location.href = "chart.html";*/
			event.preventDefault();
          });
		  
		  function plot_chart(result){
						google.charts.load('current', {packages: ['corechart', 'line']});
						google.charts.setOnLoadCallback(drawBasic);

						function drawBasic() {

							  var data = new google.visualization.DataTable();
							  data.addColumn('number', 'Curvature');
							  data.addColumn('number', 'Moment');

							  data.addRows(result[0].length);
							  for (var i=0; i<result[0].length; i++){	  
								data.setCell(i,0,result[0][i]);
								data.setCell(i,1,result[1][i]);
							}
							
						  var options = {
							legend: 'none',
							hAxis: {
							  title: 'Curvature'
							},
							vAxis: {
							  title: 'Moment'
							},
						  };

				  var chart = new google.visualization.LineChart(document.getElementById('chart_div'));

				  chart.draw(data, options);
				}
		  }
          });