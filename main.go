package main

import (
	"fmt"
	"image/color"
	"math"
	"strconv"

	"github.com/markcheno/go-quote"
	r2 "golang.org/x/exp/rand"
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/stat"
	"gonum.org/v1/gonum/stat/distmv"
	"gonum.org/v1/gonum/stat/distuv"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/plotutil"
	"gonum.org/v1/plot/vg"
)

var (
	seed          = 123457
	P     float64 = 10000.     //principle
	S01   float64 = 22.78      //GDX stock price 2018/05/14
	S02   float64 = 36.78      //XME stock price 2018/05/14
	K1    float64 = 19.8186    //S01*0.87 stike price
	K2    float64 = 31.9986    //S02*0.87 stike price
	r     float64 = 0.022      //1-Year Treasury Bill 2018/05/14
	R     float64 = 0.07 / 12. // monthly coupon rate
	path  int     = 10000      //simulation path
	N             = 240        //20 days a month, 12 months a year
	t             = 1. / 240.  //delta t
	dt            = 5          //duration (days) between issue day and settlement day
	dm            = 20         //days in a month
	month         = [12]string{"Jun", "Jul", "Aug", "Sept",
		"Oct", "Nov", "Dec", "Jan", "Feb", "March", "Apr", "May"}
)

type point struct {
	mean, std float64
}

type XYer interface {
	// Len returns the number of x, y pairs.
	Len() int

	// XY returns an x, y pair.
	XY(int) (x, y float64)
}

func main() {

	//fetch historical stock price
	GDX := Fetch_Stock("GDX")
	XME := Fetch_Stock("XME")

	//log return
	y1 := make([]float64, len(GDX)-1)
	y2 := make([]float64, len(XME)-1)
	for i := range y1 {
		y1[i] = math.Log(GDX[i+1] / GDX[i])
		y2[i] = math.Log(XME[i+1] / XME[i])
	}

	//estimate sigma
	sigma1 := stat.StdDev(y1, nil) * math.Sqrt(252)
	sigma2 := stat.StdDev(y2, nil) * math.Sqrt(252)
	fmt.Println("sigma1:", sigma1)
	fmt.Println("sigma2:", sigma2)

	//estimate correlation coefficient
	rho := stat.Correlation(y1, y2, nil)
	fmt.Println("rho:", rho)

	//scatter of log return
	nn := len(y1)
	pts := make(plotter.XYs, nn)

	for i := 0; i < nn; i++ {
		pts[i].X = y1[i]
		pts[i].Y = y2[i]
	}
	min := -0.1
	max := 0.1
	PlotScatter(pts, nn, min, max)

	// two assets option pricing random numbers
	mu := []float64{0., 0.}
	data := []float64{1, rho, rho, 1}
	cov := mat.NewSymDense(2, data)
	src := r2.New(r2.NewSource(uint64(seed)))
	normal, _ := distmv.NewNormal(mu, cov, src)

	//simulate one path and check if it meets with the rho
	e, d := AssetPrice(sigma1, sigma2, rho, normal)
	PlotGraph(e, "GDX_Sim")
	PlotGraph(d, "XME_Sim")

	//simulate 10000 paths to evaluate the product
	a := make([]float64, N+dt)    //predicted stock price for asset1
	b := make([]float64, N+dt)    //predicted stock price for asset2
	v := make([]float64, path)    //present value of all path
	v1 := make([]float64, path)   // v[i]=1 for earlier call;v[i]=0 for no earlier call
	c := make([][]float64, path)  // monthly coupon of all path
	count1 := make([]float64, 12) //number of earlier call for all paths
	count2 := 0.                  //days during 1 mo both close prices >= K1, K2
	count3 := 0.
	count4 := 0.
	sum := 0.
	sum1 := 0.
	sum2 := 0.

	for i := 0; i < path; i++ {
		a, b = AssetPrice(sigma1, sigma2, rho, normal)
		c[i] = make([]float64, 12)
		check := 0

		for j := 0; j <= 11; j++ {
			//check earlier call
			if check == 1 {
				break
			}

			count2 = 0.

			//number of days during 1 month
			//that both close prices  >= K1, K2
			for m := 1; m <= dm; m++ {
				if a[j*dm+dt+m-1] >= K1 && b[j*dm+dt+m-1] >= K2 {
					count2 = count2 + 1
				}
			}

			//monthly coupon
			if j == 0 {
				c[i][j] = P * R
			} else {
				c[i][j] = P * R * count2 / float64(dm)
			}

			//discount factor
			yy := (float64(j+1)*float64(dm) + float64(dt) - 1) / float64(N)
			disc := math.Pow((1. + r), yy)

			//present value for earlier call or non-earlier call
			if a[(j+1)*dm+dt-1] >= S01 && b[(j+1)*dm+dt-1] >= S02 { //check earlier call
				if j < 11 {
					check = 1
					count1[j] = count1[j] + 1
					v[i] = v[i] + P/disc + c[i][j]/disc
					v1[i] = 1
				} else { //no earlier call for the last month
					if a[(j+1)*dm+dt-1] < K1 || b[(j+1)*dm+dt-1] < K2 {
						if (a[(j+1)*dm+dt-1] / S01) < (b[(j+1)*dm+dt-1] / S02) {
							v[i] = v[i] + P/K1*a[(j+1)*dm+dt-1]/disc + c[i][j]/disc
						} else {
							v[i] = v[i] + P/K2*b[(j+1)*dm+dt-1]/disc + c[i][j]/disc
						}
					} else {
						v[i] = v[i] + c[i][j]/disc + P/disc
					}

				}
			} else {
				if j == 11 { //last month evaluation
					if a[(j+1)*dm+dt-1] < K1 || b[(j+1)*dm+dt-1] < K2 {
						if (a[(j+1)*dm+dt-1] / S01) < (b[(j+1)*dm+dt-1] / S02) {
							v[i] = v[i] + P/K1*a[(j+1)*dm+dt-1]/disc + c[i][j]/disc
						} else {
							v[i] = v[i] + P/K2*b[(j+1)*dm+dt-1]/disc + c[i][j]/disc
						}
					} else {
						v[i] = v[i] + c[i][j]/disc + P/disc
					}

				} else { //evaluation for every month without no ealier call
					v[i] = v[i] + c[i][j]/disc
				}
			}
		}
	}

	/*
		//test a single path
		fmt.Println(v[1])
		for j := 0; j <= 11; j++ {
			fmt.Println(c[1][j])
		}
	*/

	//(a)
	//Product Price Evaluation
	mean, std := stat.MeanStdDev(v, nil)
	std1 := stat.StdErr(std, float64(path)) //std1 = std / sqrt(path)
	results := point{mean, std1}
	fmt.Println("Product Price Evaluation(Mean,STD):", results)

	//(b)
	//print probability of earlier call for every month
	for i := 0; i < 12; i++ {
		monthindex := month[i]
		fmt.Println(monthindex, "'s probability of earlier call:", count1[i]/float64(path))
	}

	//(e)(f)
	//print probability and present value of earlier call/no earlier call
	for i := 0; i < path; i++ {
		sum = sum + v1[i]
		if v1[i] == 1 {
			sum1 = sum1 + v[i]
			count3 = count3 + 1
		} else {
			sum2 = sum2 + v[i]
			count4 = count4 + 1
		}
	}
	fmt.Println("Probability of earlier call:", sum/float64(path))
	fmt.Println("Probability of no earlier call:", 1-sum/float64(path))
	fmt.Println("Present value of earlier call:", sum1/float64(count3))
	fmt.Println("Present value of no earlier call:", sum2/float64(count4))
	fmt.Println("(b)x(e)+(1-(b))x(f)=",
		sum/float64(path)*sum1/float64(count3)+(1-sum/float64(path))*sum2/float64(count4))
	count3 = 0.
	count4 = 0.

	//print monthly coupon
	for i := 0; i < 12; i++ {
		for j := 0; j < path; j++ {
			count3 = count3 + c[j][i]
			if v[j] >= P {
				count4 = count4 + 1
			}
		}
		monthindex := month[i]
		fmt.Println(monthindex, "'s coupon:", count3/float64(path))
		count3 = 0.
	}

	//(c)
	//print probability without loss
	fmt.Println("Probability without loss:", count4/float64(path)/12.)

	//(d) VaR
	histogram(v)

	/* 錯誤的算法
	VaR95 := 1.645 * std
	VaR99 := 2.326 * std
	fmt.Println("95% chance of return:", 10000-VaR95)
	fmt.Println("99% chance of return:", 10000-VaR99)
	*/

	sum3 := make([]float64, 3000)

	for i := 0; i < 3000; i++ {
		for j := range v {
			if v[j] >= float64(i)+6000. {
				sum3[i] = sum3[i] + 1
			}
		}
	}

	for i := range sum3 {
		if sum3[i] == 9900 {
			fmt.Println("99% chance of return:", i+6000)
		}
		if sum3[i] == 9500 {
			fmt.Println("95% chance of return:", i+6000)
		}

	}

}

func Fetch_Stock(index string) []float64 {
	// fetchData
	sp, _ := quote.NewQuoteFromYahoo(index, "2017-05-15", "2018-05-14", quote.Daily, true) //AdjClose
	//fmt.Println(len(sp.Close))   // len:251  format:[]float64
	//PlotGraph(sp.Close, index)
	return sp.Close
}

func AssetPrice(sigma1, sigma2, rho float64, normal *distmv.Normal) ([]float64, []float64) {
	dd := N + dt
	s1 := make([]float64, dd) //log return  of asset1 price dist *distmv.Normal
	s2 := make([]float64, dd) //log return  of asset2 price dist *distmv.Normal
	q1 := make([]float64, dd) //asset1 price
	q2 := make([]float64, dd) //asset2 price
	s1[0] = math.Log(S01)
	s2[0] = math.Log(S02)
	q1[0] = S01
	q2[0] = S02

	z := make([][]float64, dd)

	for i := 1; i < dd; i++ {
		z[i] = make([]float64, 2)
		z[i] = normal.Rand(nil)
		s1[i] = s1[i-1] + (r-math.Pow(sigma1, 2)/2)*t + sigma1*math.Sqrt(t)*z[i][0]
		q1[i] = math.Exp(s1[i])
		s2[i] = s2[i-1] + (r-math.Pow(sigma2, 2)/2)*t + sigma2*math.Sqrt(t)*z[i][1]
		q2[i] = math.Exp(s2[i])
	}

	return q1, q2
}

func PlotGraph(a []float64, b string) {
	p, err := plot.New()
	if err != nil {
		panic(err)
	}

	p.Title.Text = "Predicted Stock Price"
	p.X.Label.Text = "Date"
	p.Y.Label.Text = "Price"
	p.Add(plotter.NewGrid())

	m := len(a)
	pts := make(plotter.XYs, m)

	for i := range pts {
		pts[i].X = float64(i)
		pts[i].Y = a[i]
	}

	// Make a line plotter with points and set its style.
	lpLine1, err := plotter.NewLine(pts)
	if err != nil {
		panic(err)
	}

	lpLine1.Color = color.RGBA{G: 255, A: 255}

	// Add the plotters to the plot
	p.Add(lpLine1)

	// Save the plot to a PNG file.
	if err := p.Save(8*vg.Inch, 6*vg.Inch, (b + ".png")); err != nil {
		panic(err)
	}

}

func histogram(x []float64) *plot.Plot {

	// convert y's type from []float64 to plotter.Values
	v := make(plotter.Values, len(x))
	for i := range v {
		v[i] = x[i]
	}
	// Make a plot and set its title.
	p, err := plot.New()
	if err != nil {
		panic(err)
	}
	p.Title.Text = "Histogram"

	// Create a histogram of our values drawn from the standard normal.
	h, err := plotter.NewHist(v, 20)
	if err != nil {
		panic(err)
	}
	h.LineStyle.Width = vg.Length(2)
	h.Color = color.RGBA{R: 255, B: 128, A: 255} // or h.LineStyle.Color = plotutil.Color(5)
	h.FillColor = plotutil.Color(6)

	// Normalize the area under the histogram to sum to one.
	h.Normalize(1)
	p.Add(h)

	// The normal distribution function
	mean := stat.Mean(v, nil)
	stdev := stat.StdDev(v, nil)
	dist := distuv.Normal{
		Mu:    mean,
		Sigma: stdev,
	}

	norm := plotter.NewFunction(dist.Prob)
	norm.Color = color.RGBA{R: 160, G: 32, B: 240, A: 255} //A:opacity
	norm.Width = vg.Points(3)
	p.Add(norm)

	// Save the plot to a PNG file.
	if err := p.Save(4*vg.Inch, 4*vg.Inch, "NormHist.png"); err != nil {
		panic(err)
	}
	return p
}

func PlotScatter(xys XYer, n int, min float64, max float64) *plot.Plot {
	// Create a new plot, set its title and  axis labels.
	p, err := plot.New()
	if err != nil {
		panic(err)
	}
	p.Title.Text = "Scatter Graph"
	str1 := strconv.Itoa(n)
	p.X.Label.Text = "n=" + str1

	// Draw a grid behind the data
	p.Add(plotter.NewGrid())

	// Make a scatter plotter and set its style.
	s, err := plotter.NewScatter(xys)
	if err != nil {
		panic(err)
	}
	s.GlyphStyle.Color = color.RGBA{R: 255, B: 128, A: 255}

	// Add the plotters to the plot, with a legend entry
	p.Add(s)
	//p.Legend.Add("scatter", s)

	//set the axis ranges
	p.X.Min = min
	p.X.Max = max

	// Save the plot to a PNG file     (6*vg.Inch)
	if err := p.Save(540, 540, "ScatterPoints.png"); err != nil {
		panic(err)
	}
	return p
}
