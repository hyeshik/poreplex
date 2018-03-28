/*
 * npinterface.cc -
 * A Python extension module to interface some functions of nanopolish.
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


#include "Python.h"

#include "nanopolish_common.h"
#include "nanopolish_squiggle_read.h"

static PyObject *ErrorObject;

PyDoc_STRVAR(npinterface_get_calibration_doc,
"get_calibration(fast5, sequence)\n\
\n\
Return the signal calibration parameters from nanopolish.");

static PyObject *
npinterface_get_calibration(PyObject *self, PyObject *args)
{
    char *fast5, *sequence;
    int res;

    if (!PyArg_ParseTuple(args, "ss:get_calibration",
                          &fast5, &sequence))
        return NULL;

    ReadDB readdb;
    std::string readid="octopus";
    std::string fast5_cpp=fast5;
    std::string sequence_cpp=sequence;

    readdb.add_signal_path(readid, fast5_cpp);
    readdb.set_intercepting_sequence(sequence_cpp);

    SquiggleRead sr(readid, readdb, SRF_LOAD_RAW_SAMPLES);

    return Py_BuildValue("ndddddd", (Py_ssize_t)sr.events[0].size(),
                         sr.scalings[0].scale, sr.scalings[0].shift,
                         sr.scalings[0].drift, sr.scalings[0].var,
                         sr.scalings[0].scale_sd, sr.scalings[0].var_sd);
}

/* List of functions defined in the module */

static PyMethodDef npinterface_methods[] = {
    {"get_calibration", npinterface_get_calibration, METH_VARARGS,
        npinterface_get_calibration_doc},
    {NULL,              NULL}           /* sentinel */
};

PyDoc_STRVAR(module_doc,
"Python interfaces to internal functions in nanopolish.");

/* Initialization function for the module (*must* be called PyInit_npinterface) */


static struct PyModuleDef npinterfacemodule = {
    PyModuleDef_HEAD_INIT,
    "npinterface",
    module_doc,
    -1,
    npinterface_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC
PyInit_npinterface(void)
{
    PyObject *m = NULL;

    /* Create the module and add the functions */
    m = PyModule_Create(&npinterfacemodule);
    if (m == NULL)
        goto fail;

    /* Add some symbolic constants to the module */
    if (ErrorObject == NULL) {
        ErrorObject = PyErr_NewException("npinterface.error", NULL, NULL);
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
