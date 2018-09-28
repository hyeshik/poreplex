/*
 * csupport.c -
 * C implementations for various poreplex functions.
 *
 *
 * Copyright (c) 2018 Hyeshik Chang
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

#define NPY_NO_DEPRECATED_API NPY_1_12_API_VERSION

#include "Python.h"
#include "numpy/arrayobject.h"
#include "event_detection.h"

static PyObject *ErrorObject;
static PyArray_Descr *event_list_descr;

static PyObject *
build_event_array(event_table *et)
{
    PyArrayObject *evarr;
    void *datap;
    npy_intp dim, i, elemsize, stride;
    event_t *evp;

    dim = et->n;

    Py_INCREF(event_list_descr);
    evarr = (PyArrayObject *)PyArray_SimpleNewFromDescr(1, &dim, event_list_descr);
    if (evarr == NULL)
        return NULL;

    /* numpy array maintains less padding than the standard C struct's.
     * Copy the data element-by-element individually. */
    elemsize = event_list_descr->elsize;
    stride = PyArray_STRIDE(evarr, 0);
    datap = PyArray_DATA(evarr);
    evp = et->event;
    for (i = 0; i < dim; i++) {
        memcpy(datap, evp++, elemsize);
        datap += stride;
    }

    return (PyObject *)evarr;
}

PyDoc_STRVAR(poreplex_detect_events_doc,
"detect_events(signal)\n\
\n\
Detect the events in the raw signal array.");

static PyObject *
poreplex_detect_events(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject *signalarg, *ret;
    PyArrayObject *signalobj;
    npy_intp nsamples;
    raw_table rt;
    event_table et;
    detector_param edparams={30, 120, 3.0f, 9.0f, 8.0f};

    static char *keywords[] = {
        "signal", "window_length1", "window_length2", "threshold1",
        "threshold2", "peak_height", NULL
    };

    assert(sizeof(Py_ssize_t) == sizeof(ssize_t));
    if (PyArg_ParseTupleAndKeywords(args, kw, "O|nnfff:detect_events", keywords,
                                    &signalarg, &edparams.window_length1,
                                    &edparams.window_length2, &edparams.threshold1,
                                    &edparams.threshold2, &edparams.peak_height) == 0)
        return NULL;

    signalobj = (PyArrayObject *)PyArray_FROM_OTF(signalarg, NPY_FLOAT,
                                 NPY_ARRAY_IN_ARRAY | NPY_ARRAY_FORCECAST);
    if (signalobj == NULL)
        return NULL;

    if (PyArray_NDIM(signalobj) != 1) {
        Py_DECREF(signalobj);
        PyErr_SetString(PyExc_ValueError, "Expects an 1-dimensional array.");
        return NULL;
    }

    nsamples = PyArray_SHAPE(signalobj)[0];

    rt.raw = PyArray_DATA(signalobj);
    rt.n = rt.end = nsamples;
    rt.start = 0;

    Py_BEGIN_ALLOW_THREADS
    et = detect_events(rt, edparams);
    Py_END_ALLOW_THREADS
    if (et.n <= 0) {
        PyErr_SetString(ErrorObject, "Event detection failed.");
        return NULL;
    }

    Py_DECREF(signalobj);

    ret = build_event_array(&et);
    if (et.event != NULL)
        free(et.event);

    return ret;
}

/* List of functions defined in the module */

static PyMethodDef poreplex_methods[] = {
    {"detect_events", (PyCFunction)poreplex_detect_events, METH_VARARGS | METH_KEYWORDS,
        poreplex_detect_events_doc},
    {NULL,              NULL}           /* sentinel */
};

PyDoc_STRVAR(module_doc,
"Supporting functions for poreplex.");

/* Initialization function for the module (*must* be called PyInit_poreplex) */

static struct PyModuleDef poreplexmodule = {
    PyModuleDef_HEAD_INIT,
    "poreplex.csupport",
    module_doc,
    -1,
    poreplex_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

static int
init_event_detector(void)
{
    PyObject *op;

#if SIZEOF_INT == 4
    op = Py_BuildValue("[(s, s), (s, s), (s, s), (s, s), (s, s), (s, s)]",
        "start", "u8", "length", "f4", "mean", "f4", "stdv", "f4",
        "pos", "i4", "state", "i4");
#elif SIZEOF_INT == 8
    /* When sizeof(int) == 8, there should be a padding for alignment between stdv and pos */
    op = Py_BuildValue("[(s, s), (s, s), (s, s), (s, s), (s, s), (s, s)]",
        "start", "u8", "length", "f4", "mean", "f4", "stdv", "f4",
        "_dummy", "i4", "pos", "i8", "state", "i8");
#else
#warning Unsupported system
#endif
    if (op == NULL)
        return -1;

    if (PyArray_DescrConverter(op, &event_list_descr) < 0) {
        Py_DECREF(op);
        return -1;
    }

    Py_DECREF(op);
    return 0;
}

PyMODINIT_FUNC
PyInit_csupport(void)
{
    PyObject *m = NULL;

    import_array();
    if (init_event_detector() < 0)
        goto fail;

    /* Create the module and add the functions */
    m = PyModule_Create(&poreplexmodule);
    if (m == NULL)
        goto fail;

    /* Add some symbolic constants to the module */
    if (ErrorObject == NULL) {
        ErrorObject = PyErr_NewException("poreplex.csupport.error", NULL, NULL);
        if (ErrorObject == NULL)
            goto fail;
    }
    Py_INCREF(ErrorObject);
    PyModule_AddObject(m, "error", ErrorObject);

    return m;
 fail:
    Py_XDECREF(m);
    return NULL;
}
