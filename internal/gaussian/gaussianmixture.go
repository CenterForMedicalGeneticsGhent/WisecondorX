package gaussian

// Copyright (c) 2014 AKUALAB INC., All rights reserved.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

import (
	"encoding/json"
	"fmt"
	"io"
	"io/ioutil"
	"math"
	"math/rand"
	"os"
	"path/filepath"
	"reflect"

	"github.com/akualab/gjoa/floatx"
	"github.com/akualab/gjoa/model"
	"github.com/akualab/gjoa/model/gaussian"
	"github.com/golang/glog"
	"github.com/gonum/floats"
)

// Model is a mixture of Gaussian distributions.
type Model struct {
	Type         string            `json:"type"`
	ModelName    string            `json:"name"`
	ModelDim     int               `json:"dim"`
	NSamples     float64           `json:"nsamples"`
	Diag         bool              `json:"diag"`
	NComponents  int               `json:"num_components"`
	PosteriorSum []float64         `json:"posterior_sum,omitempty"`
	Weights      []float64         `json:"-"`
	LogWeights   []float64         `json:"weights,omitempty"`
	Likelihood   float64           `json:"likelihood"`
	Components   []*gaussian.Model `json:"components,omitempty"`
	Iteration    int               `json:"iteration"`
	tmpProbs     []float64
}

// Option type is used to pass options to NewModel().
type Option func(*Model)

// NewModel creates a new Gaussian mixture model.
func NewModel(dim, numComponents int, options ...Option) *Model {

	gmm := &Model{
		ModelName:   "GMM", // default name
		ModelDim:    dim,
		NComponents: numComponents,
		Diag:        true,
		tmpProbs:    make([]float64, numComponents),
	}
	gmm.Type = reflect.TypeOf(*gmm).String()

	// Set options.
	for _, option := range options {
		option(gmm)
	}
	if len(gmm.PosteriorSum) == 0 {
		gmm.PosteriorSum = make([]float64, gmm.NComponents)
	}

	// Create components if not provided.
	if len(gmm.Components) == 0 {
		gmm.Components = make([]*gaussian.Model, numComponents, numComponents)
		for i := range gmm.Components {
			cname := componentName(gmm.ModelName, i, gmm.NComponents)
			gmm.Components[i] = gaussian.NewModel(gmm.ModelDim, gaussian.Name(cname))
		}
	}

	// Initialize weights.
	// Caller may pass weight, log(weights), or no weights.
	switch {

	case len(gmm.LogWeights) > 0 && len(gmm.Weights) > 0:
		glog.Fatal("options not allowed: provide only one of LogWeights or Weights")

	case len(gmm.LogWeights) == 0 && len(gmm.Weights) == 0:
		gmm.LogWeights = make([]float64, numComponents)
		logw := -math.Log(float64(gmm.NComponents))
		floatx.Apply(floatx.SetValueFunc(logw), gmm.LogWeights, nil)
		gmm.Weights = make([]float64, gmm.NComponents)
		floatx.Exp(gmm.Weights, gmm.LogWeights)
		glog.Infof("init weights with equal values: %.6f", gmm.Weights[0])

	case len(gmm.LogWeights) > 0:
		gmm.Weights = make([]float64, gmm.NComponents)
		floatx.Exp(gmm.Weights, gmm.LogWeights)

	case len(gmm.Weights) > 0:
		gmm.LogWeights = make([]float64, numComponents)
		floatx.Log(gmm.LogWeights, gmm.Weights)
	}
	return gmm
}

// Computes log prob for mixture.// SIDE EFFECT => returns logProb of
// Gaussian comp + logWeight in matrix pointed by func arg probs.
func (gmm *Model) logProbInternal(obs, probs []float64) float64 {

	var max = -math.MaxFloat64

	/* Compute log probabilities for this observation. */
	for i, c := range gmm.Components {
		o := model.F64ToObs(obs, "")
		v1 := c.LogProb(o)
		v2 := gmm.LogWeights[i]
		v := v1 + v2

		if probs != nil {
			probs[i] = v
		}
		if v > max {
			max = v
		}
	}

	// To simplify computation, use the max prob in the denominator instead
	// of the sum.
	return max
}

// LogProb returns log probability for observation.
func (gmm *Model) LogProb(obs model.Obs) float64 {

	o := obs.Value().([]float64)
	return gmm.logProbInternal(o, nil)
}

// Returns the probability.
func (gmm *Model) prob(obs []float64) float64 {
	return math.Exp(gmm.LogProb(model.F64ToObs(obs, "")))
}

// Predict returns a hypothesis given the observation.
func (gmm *Model) Predict(x model.Observer) ([]model.Labeler, error) {

	glog.Fatal("Predict method not implemented.")
	return nil, nil
}

/*
  The posterior prob for each mixture component. We approximate the sum using max.

                 p(o|c(i)) p(c(i))
   p(c(i)|o) ~ ---------------------
                max{p(o|c(i)) p(c(i))}

*/

// Estimate computes model parameters using sufficient statistics.
func (gmm *Model) UpdateOne(o model.Obs, w float64) {

	obs, _, _ := model.ObsToF64(o)
	maxProb := gmm.logProbInternal(obs, gmm.tmpProbs)
	gmm.Likelihood += maxProb
	floatx.Apply(floatx.AddScalarFunc(-maxProb+math.Log(w)), gmm.tmpProbs, nil)

	// Compute posterior probabilities.
	floatx.Exp(gmm.tmpProbs, gmm.tmpProbs)

	// Update posterior sum, needed to compute mixture weights.
	floats.Add(gmm.PosteriorSum, gmm.tmpProbs)

	// Update Gaussian components.
	for i, c := range gmm.Components {
		c.UpdateOne(o, gmm.tmpProbs[i])
	}

	// Count number of observations.
	gmm.NSamples += w
}

// Update updates sufficient statistics using observations.
func (gmm *Model) Update(x model.Observer, w func(model.Obs) float64) error {
	c, e := x.ObsChan()
	if e != nil {
		return e
	}
	for v := range c {
		gmm.UpdateOne(v, w(v))
	}
	return nil
}

// Estimate computes model parameters using sufficient statistics.
func (gmm *Model) Estimate() error {

	// Estimate mixture weights.
	floatx.Apply(floatx.ScaleFunc(1.0/gmm.NSamples), gmm.PosteriorSum, gmm.Weights)
	floatx.Log(gmm.LogWeights, gmm.Weights)

	// Estimate component density.
	for _, c := range gmm.Components {
		err := c.Estimate()
		if err != nil {
			return err
		}
	}
	gmm.Iteration++

	return nil
}

// Clear resets sufficient statistics.
func (gmm *Model) Clear() {

	for _, c := range gmm.Components {
		c.Clear()
	}
	floatx.Apply(floatx.SetValueFunc(0), gmm.PosteriorSum, nil)
	gmm.NSamples = 0
	gmm.Likelihood = 0
}

// Sample returns a GMM sample.
func (gmm *Model) Sample(r *rand.Rand) model.Obs {
	// Choose a component using weights
	comp := model.RandIntFromDist(gmm.Weights, r)
	// Get a random vector from that component
	return gmm.Components[comp].Sample(r)
}

// SampleChan returns a channel with samples generated by the GMM model.
func (gmm *Model) SampleChan(r *rand.Rand, size int) <-chan model.Obs {

	if len(gmm.Weights) == 0 {
		glog.Fatal("Parameter Weights is missing.")
	}
	c := make(chan model.Obs, 1000)
	go func() {
		for i := 0; i < size; i++ {
			c <- gmm.Sample(r)
		}
		close(c)
	}()
	return c
}

// Returns a random vector using the mean and sd vectors.
func RandomVector(mean, sd []float64, r *rand.Rand) []float64 {

	nrows := len(mean)
	if !floats.EqualLengths(mean, sd) {
		panic(floatx.ErrLength)
	}

	vec := make([]float64, nrows)
	for i := 0; i < nrows; i++ {
		v := r.NormFloat64()*sd[i] + mean[i]
		vec[i] = v
	}
	return vec
}

// RandomModel generates a random Gaussian mixture model using mean and variance vectors as seed.
// Use this function to initialize the GMM before training. The mean and sd
// vector can be estimated from the data set using a Gaussian model.
func RandomModel(mean, sd []float64, numComponents int,
	name string, seed int64) *Model {

	n := len(mean)
	if !floats.EqualLengths(mean, sd) {
		panic(floatx.ErrLength)
	}
	cs := make([]*gaussian.Model, n, n)
	r := rand.New(rand.NewSource(seed))
	for i := 0; i < n; i++ {
		rv := RandomVector(mean, sd, r)
		cs[i] = gaussian.NewModel(n, gaussian.Mean(rv), gaussian.StdDev(sd))
	}
	gmm := NewModel(n, numComponents, Name(name), Components(cs))
	return gmm
}

func componentName(name string, n, numComponents int) string {

	max := numComponents - 1
	switch {
	case max < 10:
		return fmt.Sprintf("%s-%d", name, n)

	case max < 100:
		return fmt.Sprintf("%s-%02d", name, n)

	case max < 1000:
		return fmt.Sprintf("%s-%03d", name, n)

	case max < 10000:
		return fmt.Sprintf("%s-%04d", name, n)

	default:
		return fmt.Sprintf("%s-%d", name, n)
	}
}

// Dim is the dimensionality of the observation vector.
func (gmm *Model) Dim() int { return gmm.ModelDim }

// Export struct.
/*type GMMValues struct {
	Name          string
	Type          string
	NumElements   int
	NumSamples    float64
	Diagonal      bool
	NumComponents int
	PosteriorSum  []float64 `json:",omitempty"`
	Weights       []float64
	Likelihood    float64
	Components    []interface{}
}

func (gmm *GMM) Values() interface{} {

	values := &GMMValues{
		Name:          gmm.name,
		Type:          gmm.Type(),
		NumElements:   gmm.numElements,
		NumSamples:    gmm.numSamples,
		Diagonal:      gmm.diagonal,
		NumComponents: gmm.numComponents,
		PosteriorSum:  gmm.PosteriorSum,
		Weights:       gmm.Weights,
		Likelihood:    gmm.totalLikelihood,
		Components:    make([]interface{}, gmm.numComponents),
	}

	for k, v := range gmm.components {
		values.Components[k] = v.Values()
	}

	return values
}
*/
/*
func (gmm *GMM) New(values interface{}) (model.Modeler, error) {

	v := values.(*GMMValues)

	ng, e := NewGaussianMixture(v.NumElements, v.NumComponents, true, v.Diagonal, v.Name)
	if e != nil {
		return nil, e
	}

	gaussian := &Gaussian{}
	for k, gv := range v.Components {

		// Create a Gaussian component using the Gaussian values.
		g, e := gaussian.New(gv)
		if e != nil {
			return nil, e
		}
		ng.components[k] = g.(*Gaussian)
	}
	ng.Weights = v.Weights
	ng.numSamples = v.NumSamples
	ng.totalLikelihood = v.Likelihood

	if len(v.PosteriorSum) > 0 {
		ng.PosteriorSum = v.PosteriorSum
	}

	return ng, nil
}
*/

// Name returns the name of the model.
func (gmm *Model) Name() string {
	return gmm.ModelName
}

// SetName sets a name for the model.
func (gmm *Model) setName(name string) {
	gmm.ModelName = name
}

// Options

// Name is an option to set the model name.
func Name(name string) Option {
	return func(gmm *Model) {
		glog.Infof("set model name: %s", gmm.Name)
		gmm.setName(name)
	}
}

// Components sets the mixture components for the model.
func Components(cs []*gaussian.Model) Option {
	return func(gmm *Model) {
		gmm.Components = cs
		glog.Infof("set %d mixture components.", len(gmm.Components))
	}
}

// Weights sets the mixture weights for the model.
func Weights(w []float64) Option {
	return func(gmm *Model) { gmm.Weights = w }
}

// LogWeights sets the mixture weights for the model
// using log(w) as the argument.
func LogWeights(logw []float64) Option {
	return func(gmm *Model) { gmm.LogWeights = logw }
}

// Clone create a clone of src.
func Clone(src *Model) Option {
	return func(m *Model) {
		m.NSamples = src.NSamples
		m.Diag = src.Diag
		m.PosteriorSum = src.PosteriorSum
		m.Iteration = src.Iteration
	}
}

// IO

// Read unmarshals json data from an io.Reader into a model struct.
func Read(r io.Reader) (*Model, error) {

	b, err := ioutil.ReadAll(r)
	if err != nil {
		return nil, err
	}

	// Get a Model object.
	m := &Model{}
	e := json.Unmarshal(b, m)

	if e != nil {
		return nil, e
	}
	m = NewModel(m.ModelDim, m.NComponents, Clone(m), LogWeights(m.LogWeights),
		Components(m.Components), Name(m.ModelName))
	return m, nil
}

// ReadFile unmarshals json data from a file into a model struct.
func ReadFile(fn string) (*Model, error) {

	f, err := os.Open(fn)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	glog.Infof("Reading model from file %s.", fn)
	return Read(f)
}

// Write writes the model to an io.Writer.
func (gmm *Model) Write(w io.Writer) error {

	b, err := json.Marshal(gmm)
	if err != nil {
		return err
	}
	_, e := w.Write(b)
	return e
}

// WriteFile writes the model to file.
func (gmm *Model) WriteFile(fn string) error {

	e := os.MkdirAll(filepath.Dir(fn), 0755)
	if e != nil {
		return e
	}
	f, err := os.Create(fn)
	if err != nil {
		return err
	}
	defer f.Close()

	ee := gmm.Write(f)
	if ee != nil {
		return ee
	}

	glog.Infof("Wrote model \"%s\" to file %s.", gmm.Name(), fn)
	return nil
}
