from shogun.Features import *
from shogun.Kernel import *
from shogun.Classifier import *
from shogun.Evaluation import *
from shogun.Library import *
from numpy.random import rand
import inputParser
import parse
import numpy
import sys

label_test = [] #Global list to handle test data.

def initClassifiers(inputDir):
    """ This method is the starting point for SVM implementation """
    result = inputParser.getInput(inputDir)
    label = result['label']

    print 'Pruning the dataset..'    

    totalCount = 70000  # total no of samples to be considered.
    trainCount = (totalCount * 80)/100 # 80% of total count used as training data.
    testCount = totalCount - trainCount # rest 20% is used as testing data.

    print 'Total no of samples considered: %d' % totalCount
    print 'No of samples trained on: %d' % trainCount
    print 'No of testing samples: %d' % (testCount)

    # Inputs - sliced according to size.
    train_input = (result['inputs'])[:trainCount]
    label_train_H = numpy.asarray((result['label_H'])[:trainCount])
    label_train_E = numpy.asarray((result['label_E'])[:trainCount])
    label_train_C = numpy.asarray((result['label_C'])[:trainCount])
    
    test_input = (result['inputs'])[trainCount:totalCount]
    # Labels for each individual classes.
    label_test_H = (result['label_H'])[trainCount:totalCount]
    label_test_E = (result['label_E'])[trainCount:totalCount]
    label_test_C = (result['label_C'])[trainCount:totalCount]
    
    # Label list. - contains labels for all the classes.
    label_test.extend((result['label'])[trainCount:totalCount])

    #Orthogonalized input.
    train_input_feats = numpy.transpose(numpy.asarray((result['featMatrix'])[:trainCount]))
    test_input_feats = numpy.transpose(numpy.asarray((result['featMatrix'])[trainCount:totalCount]))
    label_train = numpy.asarray((result['label'])[:trainCount])
 
    # parameters for SVMLight implementation.
#    parameter = [train_input,test_input,label_train_H,label_train_E,label_train_C, label_test_H,label_test_E,label_test_C, 5, 1]
                 
#    classifier_svmlight(*parameter)
    
    # Invoking Experimental SVMLib implementation.
    out_new = bayesianClassifier(train_input_feats,test_input_feats, label_train, label_test)
    print out_new
    output = '--'
    for i in range(0, len(out_new)):
        if out_new[i] == 1 :
            output += 'H'
        elif out_new[i] == 2 :
            output += 'E'
        else :
            output += 'C'
    print output + '--'


def classifier_svmlight(fm_train_input,fm_test_inp,
                               label_train_H, label_train_E,
                               label_train_C, label_test_H, 
                               label_test_E, label_test_C, 
                               degree=3, C=10):
    """ Classifier using SVMLight library and string kernel """

    # Defining features using shogun library.
    feats_train=StringCharFeatures(PROTEIN)
    feats_train.set_features(fm_train_input)
   
    feats_test=StringCharFeatures(PROTEIN)
    feats_test.set_features(fm_test_inp)
    
    # Defining labels for each class.
    labels_H = BinaryLabels(label_train_H)
    labels_E = BinaryLabels(label_train_E)
    labels_C = BinaryLabels(label_train_C)    
    
    kernel = WeightedDegreePositionStringKernel(feats_train, feats_train, degree)

    # Instantiating different SVM's for classes using same kernel.
    svm_H = SVMLight(C, kernel, labels_H)
    svm_E = SVMLight(C, kernel, labels_E)
    svm_C = SVMLight(C, kernel, labels_C)

    print 'Training on H..'
    svm_H.train()
    print 'Training on E..'
    svm_E.train()
    print 'Training on C..'
    svm_C.train()

    # Combining the three SVM's to achieve multiclass SVM using oneVsRest.
    print 'Running svm on test data..'
    out_H = svm_H.apply(feats_test).get_labels()

    label_test_E = []
    fm_test_E = []

    for i in range(0, len(out_H)):
        if label_test[i] == 1 and out_H[i] == 1:
            continue;
        label_test_E.append(label_test[i])
        fm_test_E.append(fm_test_inp[i])

    feats_test_E=StringCharFeatures(PROTEIN)
    feats_test_E.set_features(fm_test_E)

    out_E = svm_E.apply(feats_test_E).get_labels()

    label_test_C = []
    fm_test_C = []

    for i in range(0, len(out_E)):
        if label_test_E[i] == 2 and out_E[i] == 1:
            continue;
        label_test_C.append(label_test_E[i])
        fm_test_C.append(fm_test_E[i])

    feats_test_C=StringCharFeatures(PROTEIN)
    feats_test_C.set_features(fm_test_C)

    out_C = svm_C.apply(feats_test_C).get_labels()

    # Measuring the accuracy using Q3.
    count_H = 0.0
    count_miss_H = 0.0
    count_E = 0.0
    count_miss_E = 0.0
    count_C = 0.0
    count_miss_C = 0.0
    correctCount = 0.0

    print 'Analysing the results.. '

    # Stats from actual test data. 
    countH_ori = label_test.count(1)
    countE_ori = label_test.count(2)
    countC_ori = label_test.count(3)

    # Stats of correctly predicted in class H.
    for i in range(0, len(out_H)):
        if label_test[i] == 1 and out_H[i] == 1:
            count_H += 1

    # Stats of correctly predicted in class E.
    for i in range(0, len(out_E)):
        if label_test_E[i] == 2 and out_E[i] == 1:
            count_E += 1

    # Stats of correctly predicted in class C.
    for i in range(0, len(out_C)):
        if label_test_C[i] == 3 and out_C[i] == 1:
            count_C += 1


    print 'Total predicted as H: %d' % count_H
    print 'Total predicted as E: %d' % count_E
    print 'Total predicted as C: %d' % count_C

    # Calculating the Q3 Score.
    q3 = (1.0 * (count_H + count_E + count_C))/(countH_ori + countE_ori + countC_ori)
    q3 = q3 * 100
    print 'Q3: %.4f' % q3



def bayesianClassifier(inputSeq_train, inputSeq_test, labels_train, labels_test ):
    """ An experimental implementation using orthogonalization of the input.
        Gaussian Naive Bayesian Classifier is used.
    """
    print 'Gaussian Naive Bayesian Classifier.'    

    feats_train = RealFeatures(inputSeq_train)
    feats_test = RealFeatures(inputSeq_test)

    labels = MulticlassLabels(labels_train)

    print 'Initializing kernel..'
    classifier = GaussianNaiveBayes(feats_train, labels)

    print 'Training svm..'
    classifier_train = classifier.train()

    print 'Running test data..'
    label_pred = classifier.apply(feats_test)
    out = label_pred.get_labels()

    # Evaluating the accuracy.
    if label_test is not None:
        labels_test_set = MulticlassLabels(numpy.asarray(label_test))
        evaluator = MulticlassAccuracy()
        acc = evaluator.evaluate(label_pred, labels_test_set)
    #        print 'Accuracy : %.4f' % (acc * 100)
    
    # Giving the new input.
    inp = parse.run("513_distribute/7cata.all")
    inpMatrix = inputParser.createFeatureMatrix(inp)
    inp_feats = RealFeatures(numpy.transpose(numpy.asarray(inpMatrix)))
    
    # Running svm on the input.
    label_out = classifier.apply(inp_feats)
    out_new = label_out.get_labels()
    return out_new


if __name__=='__main__':
    print('Initializing Classifiers..')
    initClassifiers(sys.argv[1])
