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

typedef struct {
    PyObject_HEAD
    SquiggleRead *sr;
    ReadDB *readdb;
    PyObject *attrdict;
} SquiggleReadObject;

#define SquiggleReadObject_Check(v)     (Py_TYPE(v) == &SquiggleRead_Type)

/* SquiggleRead methods */

static void
SquiggleRead_dealloc(SquiggleReadObject *self)
{
    if (self->sr != NULL)
        delete self->sr;
    self->sr = NULL;

    if (self->readdb != NULL)
        delete self->readdb;
    self->readdb = NULL;

    Py_XDECREF(self->attrdict);
    self->attrdict = NULL;

    PyObject_Del(self);
}

static PyObject *
SquiggleRead_get_events(SquiggleReadObject *self, PyObject *args)
{
    double ev_start, ev_end;
    size_t i, totalevents;
    std::vector<SquiggleEvent> *events;
    PyObject *evlist;

    if (!PyArg_ParseTuple(args, "dd:get_events", &ev_start, &ev_end))
        return NULL;

    evlist = PyList_New(0);
    if (evlist == NULL)
        return NULL;

    events = &self->sr->events[0];
    totalevents = events->size();
    for (i = 0; i < totalevents; i++) {
        SquiggleEvent ev=events->at(i);
        PyObject *tup;

        if (ev.start_time + ev.duration <= ev_start || ev.start_time >= ev_end)
            continue;

        tup = Py_BuildValue("ffdff", ev.mean, ev.stdv, ev.start_time,
                            ev.duration, ev.log_stdv);
        if (tup == NULL ||
                PyList_Append(evlist, tup) != 0) {
            Py_XDECREF(tup);
            Py_DECREF(evlist);
            return NULL;
        }
        Py_DECREF(tup);
    }

    return evlist;
}

static PyObject *
SquiggleRead_get_read_sequence(SquiggleReadObject *self, void *closure)
{
    return PyUnicode_FromString(self->sr->read_sequence.c_str());
}

static PyMethodDef SquiggleRead_methods[] = {
    {"get_events",              (PyCFunction)SquiggleRead_get_events,
                METH_VARARGS,
                PyDoc_STR("get_events() -> list of (mean, stdv, start_time, duration, log_stdv)")},
    {NULL,                      NULL}
};

static PyObject *
SquiggleRead_get_event_count(SquiggleReadObject *self, void *closure)
{
    return PyLong_FromLong(self->sr->events[0].size());
}

static PyObject *
SquiggleRead_get_scaling_params(SquiggleReadObject *self, void *closure)
{
    return Py_BuildValue("dddddd",
                         self->sr->scalings[0].scale,
                         self->sr->scalings[0].shift,
                         self->sr->scalings[0].drift,
                         self->sr->scalings[0].var,
                         self->sr->scalings[0].scale_sd,
                         self->sr->scalings[0].var_sd);
}

static PyObject *
SquiggleRead_get_sample_rate(SquiggleReadObject *self, void *closure)
{
    return PyFloat_FromDouble(self->sr->sample_rate);
}

static PyObject *
SquiggleRead_get_sample_start_time(SquiggleReadObject *self, void *closure)
{
    return PyLong_FromLong(self->sr->sample_start_time);
}

static PyObject *
SquiggleRead_get_base_to_event_map(SquiggleReadObject *self, void *closure)
{
    PyObject *ret;
    size_t i;

    ret = PyList_New(self->sr->base_to_event_map.size());
    if (ret == NULL)
        return NULL;

    for (i = 0; i < self->sr->base_to_event_map.size(); i++) {
        EventRangeForBase *evmap=&self->sr->base_to_event_map.at(i);

        PyObject *basestart, *basestop;
        basestart = basestop = NULL;
        if ((basestart = PyLong_FromLong(evmap->indices[0].start)) == NULL ||
            (basestop = PyLong_FromLong(evmap->indices[0].stop)) == NULL) {
            Py_XDECREF(basestart); Py_XDECREF(basestop);
            Py_DECREF(ret);
            return NULL;
        }

        PyObject *idxmap=PyTuple_New(2);
        if (idxmap == NULL) {
            Py_DECREF(basestart); Py_DECREF(basestop);
            Py_DECREF(ret);
            return NULL;
        }

        PyTuple_SET_ITEM(idxmap, 0, basestart);
        PyTuple_SET_ITEM(idxmap, 1, basestop);
        PyList_SET_ITEM(ret, i, idxmap);
    }
    return ret;
}

static PyObject *
SquiggleRead_getattro(SquiggleReadObject *self, PyObject *name)
{
    return PyObject_GenericGetAttr((PyObject *)self, name);
}

static PyGetSetDef SquiggleRead_getsetlist[] = {
    {"event_count",
     (getter)SquiggleRead_get_event_count, NULL,
     NULL, NULL},
    {"read_sequence",
     (getter)SquiggleRead_get_read_sequence, NULL,
     NULL, NULL},
    {"scaling_params",
     (getter)SquiggleRead_get_scaling_params, NULL,
     NULL, NULL},
    {"sample_rate",
     (getter)SquiggleRead_get_sample_rate, NULL,
     NULL, NULL},
    {"sample_start_time",
     (getter)SquiggleRead_get_sample_start_time, NULL,
     NULL, NULL},
    {"base_to_event_map",
     (getter)SquiggleRead_get_base_to_event_map, NULL,
     NULL, NULL},
    {NULL} /* Sentinel */
};

static PyTypeObject SquiggleRead_Type = {
    /* The ob_type field must be initialized in the module init function
     * to be portable to Windows without using C++. */
    PyVarObject_HEAD_INIT(NULL, 0)
    "npinterface.SquiggleRead", /*tp_name*/
    sizeof(SquiggleReadObject), /*tp_basicsize*/
    0,                          /*tp_itemsize*/
    /* methods */
    (destructor)SquiggleRead_dealloc, /*tp_dealloc*/
    0,                          /*tp_print*/
    (getattrfunc)0,             /*tp_getattr*/
    (setattrfunc)0,             /*tp_setattr*/
    0,                          /*tp_reserved*/
    0,                          /*tp_repr*/
    0,                          /*tp_as_number*/
    0,                          /*tp_as_sequence*/
    0,                          /*tp_as_mapping*/
    0,                          /*tp_hash*/
    0,                          /*tp_call*/
    0,                          /*tp_str*/
    (getattrofunc)SquiggleRead_getattro, /*tp_getattro*/
    0,                          /*tp_setattro*/
    0,                          /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,         /*tp_flags*/
    0,                          /*tp_doc*/
    0,                          /*tp_traverse*/
    0,                          /*tp_clear*/
    0,                          /*tp_richcompare*/
    0,                          /*tp_weaklistoffset*/
    0,                          /*tp_iter*/
    0,                          /*tp_iternext*/
    SquiggleRead_methods,       /*tp_methods*/
    0,                          /*tp_members*/
    SquiggleRead_getsetlist,    /*tp_getset*/
    0,                          /*tp_base*/
    0,                          /*tp_dict*/
    0,                          /*tp_descr_get*/
    0,                          /*tp_descr_set*/
    0,                          /*tp_dictoffset*/
    0,                          /*tp_init*/
    0,                          /*tp_alloc*/
    0,                          /*tp_new*/
    0,                          /*tp_free*/
    0,                          /*tp_is_gc*/
};
/* --------------------------------------------------------------------- */

/* List of functions defined in the module */

static PyObject *
newSquiggleRead(PyObject *_, PyObject *args)
{
    SquiggleReadObject *self;
    char *fast5;

    if (!PyArg_ParseTuple(args, "s:SquiggleRead", &fast5))
        return NULL;

    self = PyObject_New(SquiggleReadObject, &SquiggleRead_Type);
    if (self == NULL)
        return NULL;

    self->attrdict = PyDict_New();
    if (self->attrdict == NULL) {
        Py_DECREF(self);
        return NULL;
    }

    self->readdb = new ReadDB();
    if (self->readdb == NULL) {
        Py_DECREF(self);
        return NULL;
    }

    std::string readid="octopus";
    std::string fast5_cpp=fast5;
    self->readdb->add_signal_path(readid, fast5_cpp);

    self->sr = new SquiggleRead(readid, *self->readdb, SRF_LOAD_RAW_SAMPLES);
    if (self->sr == NULL) {
        Py_DECREF(self);
        return NULL;
    }

    return (PyObject *)self;
}

static PyMethodDef npinterface_methods[] = {
    {"open_squiggle_read",  (PyCFunction)newSquiggleRead,    METH_VARARGS,
            PyDoc_STR("open_squiggle_read(filename) -> SquiggleRead")},
    {NULL,              NULL}           /* sentinel */
};


static int
npinterface_exec(PyObject *m)
{
    /* Finalize the type object including setting type of the new type
     * object; doing it here is required for portability, too. */
    if (PyType_Ready(&SquiggleRead_Type) < 0)
        goto fail;

    /* Add some symbolic constants to the module */
    if (ErrorObject == NULL) {
        ErrorObject = PyErr_NewException("npinterface.error", NULL, NULL);
        if (ErrorObject == NULL)
            goto fail;
    }
    Py_INCREF(ErrorObject);
    PyModule_AddObject(m, "error", ErrorObject);

    return 0;
 fail:
    Py_XDECREF(m);
    return -1;
}

/* Initialization function for the module (*must* be called PyInit_npinterface) */

static struct PyModuleDef_Slot npinterface_slots[] = {
    {Py_mod_exec, (void *)npinterface_exec},
    {0, NULL},
};

PyDoc_STRVAR(module_doc,
"Provides nanopolish internal objects.");

static struct PyModuleDef npinterfacemodule = {
    PyModuleDef_HEAD_INIT,
    "npinterface",
    module_doc,
    0,
    npinterface_methods,
    npinterface_slots,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC
PyInit_npinterface(void)
{
    return PyModuleDef_Init(&npinterfacemodule);
}
