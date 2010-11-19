/****************************************************************************
** Meta object code from reading C++ file 'mepp_component_MSDM_plugin.hxx'
**
** Created: Mon 15. Nov 14:13:54 2010
**      by: The Qt Meta Object Compiler version 62 (Qt 4.6.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "mepp_component_MSDM_plugin.hxx"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'mepp_component_MSDM_plugin.hxx' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.6.3. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_mepp_component_MSDM_plugin[] = {

 // content:
       4,       // revision
       0,       // classname
       0,    0, // classinfo
       2,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      28,   27,   27,   27, 0x0a,
      47,   27,   27,   27, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_mepp_component_MSDM_plugin[] = {
    "mepp_component_MSDM_plugin\0\0"
    "MSDM_computation()\0DistanceToColorMap()\0"
};

const QMetaObject mepp_component_MSDM_plugin::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_mepp_component_MSDM_plugin,
      qt_meta_data_mepp_component_MSDM_plugin, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &mepp_component_MSDM_plugin::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *mepp_component_MSDM_plugin::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *mepp_component_MSDM_plugin::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_mepp_component_MSDM_plugin))
        return static_cast<void*>(const_cast< mepp_component_MSDM_plugin*>(this));
    if (!strcmp(_clname, "mepp_component_plugin_interface"))
        return static_cast< mepp_component_plugin_interface*>(const_cast< mepp_component_MSDM_plugin*>(this));
    if (!strcmp(_clname, "fr.liris.MEPP.PluginInterface/1.0"))
        return static_cast< mepp_component_plugin_interface*>(const_cast< mepp_component_MSDM_plugin*>(this));
    return QObject::qt_metacast(_clname);
}

int mepp_component_MSDM_plugin::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: MSDM_computation(); break;
        case 1: DistanceToColorMap(); break;
        default: ;
        }
        _id -= 2;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
