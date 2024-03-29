// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/UI/CQProgressItemText.h,v $
//   $Revision: 1.8 $
//   $Name: Build-33 $
//   $Author: gauges $
//   $Date: 2009/02/18 20:47:31 $
// End CVS Header

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

#ifndef CQPROGRESSITEMTEXT_H
#define CQPROGRESSITEMTEXT_H

#include <qvariant.h>
//Added by qt3to4:
#include <QPixmap>
#include <Q3HBoxLayout>
#include <QLabel>
#include "CQProgressItem.h"

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <Qt3Support/Q3ButtonGroup>
#include <Qt3Support/Q3HBoxLayout>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QSpacerItem>
#include <QtGui/QWidget>
#include "CQProgressItem.h"
#include "utilities/CProcessReport.h"
#include "utilities/CVector.h"

class Ui_CQProgressItemText
  {
  public:
    QWidget *mLabel;
    Q3HBoxLayout *hboxLayout;
    QLabel *mItemName;
    QSpacerItem *mSpacer;
    QLineEdit *mValue;

    void setupUi(CQProgressItem *CQProgressItemText)
    {
      if (CQProgressItemText->objectName().isEmpty())
        CQProgressItemText->setObjectName(QString::fromUtf8("CQProgressItemText"));
      CQProgressItemText->resize(404, 33);
      QSizePolicy sizePolicy(static_cast<QSizePolicy::Policy>(3), static_cast<QSizePolicy::Policy>(0));
      sizePolicy.setHorizontalStretch(0);
      sizePolicy.setVerticalStretch(0);
      sizePolicy.setHeightForWidth(CQProgressItemText->sizePolicy().hasHeightForWidth());
      CQProgressItemText->setSizePolicy(sizePolicy);
      CQProgressItemText->setMinimumSize(QSize(200, 33));
      mLabel = new QWidget(CQProgressItemText);
      mLabel->setObjectName(QString::fromUtf8("mLabel"));
      mLabel->setGeometry(QRect(10, 10, 330, 22));
      hboxLayout = new Q3HBoxLayout(mLabel);
      hboxLayout->setSpacing(6);
      hboxLayout->setMargin(1);
      hboxLayout->setObjectName(QString::fromUtf8("hboxLayout"));
      hboxLayout->setContentsMargins(0, 0, 0, 0);
      mItemName = new QLabel(mLabel);
      mItemName->setObjectName(QString::fromUtf8("mItemName"));
      QSizePolicy sizePolicy1(static_cast<QSizePolicy::Policy>(5), static_cast<QSizePolicy::Policy>(3));
      sizePolicy1.setHorizontalStretch(0);
      sizePolicy1.setVerticalStretch(0);
      sizePolicy1.setHeightForWidth(mItemName->sizePolicy().hasHeightForWidth());
      mItemName->setSizePolicy(sizePolicy1);
      mItemName->setMinimumSize(QSize(120, 20));
      mItemName->setWordWrap(false);

      hboxLayout->addWidget(mItemName);

      mSpacer = new QSpacerItem(16, 16, QSizePolicy::Maximum, QSizePolicy::Minimum);

      hboxLayout->addItem(mSpacer);

      mValue = new QLineEdit(mLabel);
      mValue->setObjectName(QString::fromUtf8("mValue"));
      QSizePolicy sizePolicy2(static_cast<QSizePolicy::Policy>(3), static_cast<QSizePolicy::Policy>(5));
      sizePolicy2.setHorizontalStretch(0);
      sizePolicy2.setVerticalStretch(0);
      sizePolicy2.setHeightForWidth(mValue->sizePolicy().hasHeightForWidth());
      mValue->setSizePolicy(sizePolicy2);
      mValue->setMinimumSize(QSize(200, 0));

      hboxLayout->addWidget(mValue);

      retranslateUi(CQProgressItemText);

      QMetaObject::connectSlotsByName(CQProgressItemText);
    } // setupUi

    void retranslateUi(CQProgressItem *CQProgressItemText)
    {
      CQProgressItemText->setCaption(QApplication::translate("CQProgressItemText", "Progress Text", 0, QApplication::UnicodeUTF8));
      mItemName->setText(QApplication::translate("CQProgressItemText", "Item Name", 0, QApplication::UnicodeUTF8));
      Q_UNUSED(CQProgressItemText);
    } // retranslateUi

  protected:
    enum IconID
    {
      image0_ID,
      unknown_ID
    };
    static QPixmap qt_get_icon(IconID id)
    {
      static const unsigned char image0_data[] =
        {
          0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a, 0x00, 0x00, 0x00, 0x0d,
          0x49, 0x48, 0x44, 0x52, 0x00, 0x00, 0x00, 0x16, 0x00, 0x00, 0x00, 0x16,
          0x08, 0x06, 0x00, 0x00, 0x00, 0xc4, 0xb4, 0x6c, 0x3b, 0x00, 0x00, 0x03,
          0xb1, 0x49, 0x44, 0x41, 0x54, 0x78, 0x9c, 0xad, 0x94, 0x51, 0x4c, 0x5b,
          0x55, 0x18, 0xc7, 0x7f, 0xe7, 0xdc, 0x4b, 0x7b, 0x4b, 0x61, 0x50, 0xbb,
          0x96, 0x32, 0x44, 0x18, 0xca, 0x32, 0x35, 0x8b, 0xee, 0x61, 0x92, 0x60,
          0x9c, 0x51, 0xd8, 0x83, 0x89, 0x2c, 0xe0, 0x83, 0xf1, 0x71, 0x8b, 0x3e,
          0xbb, 0x18, 0x5f, 0x8d, 0xc9, 0x1e, 0x97, 0x2c, 0xf3, 0x9d, 0x2d, 0x6a,
          0x78, 0xd0, 0xb0, 0x27, 0xb3, 0xcd, 0x07, 0xd9, 0xe6, 0x8c, 0x81, 0xc6,
          0x25, 0xa6, 0xc1, 0x39, 0x40, 0x8a, 0x43, 0x84, 0xf4, 0x16, 0x10, 0x0a,
          0xed, 0x6d, 0x4b, 0x7b, 0xdb, 0x7b, 0x8e, 0x0f, 0xc0, 0x92, 0x1a, 0x70,
          0xc6, 0xed, 0x7b, 0x3a, 0xe7, 0xe4, 0x3b, 0xbf, 0xf3, 0xff, 0xfe, 0xdf,
          0xc9, 0x27, 0xe2, 0xf1, 0x38, 0xbb, 0xd1, 0xdb, 0xdb, 0xab, 0x79, 0x02,
          0x11, 0x8f, 0xc7, 0x85, 0xd8, 0x05, 0x8f, 0x8d, 0x8d, 0xe9, 0xae, 0xae,
          0x2e, 0x72, 0xb9, 0x1c, 0xb6, 0x6d, 0xe3, 0x38, 0xce, 0x7f, 0x82, 0xe4,
          0x72, 0x39, 0x66, 0x67, 0x67, 0x49, 0x24, 0x12, 0xb8, 0xae, 0xfb, 0xf0,
          0xdc, 0xdc, 0x55, 0x3a, 0x32, 0x32, 0x42, 0xf4, 0x50, 0x1d, 0x91, 0x8e,
          0x2d, 0x3a, 0x8f, 0x15, 0xa9, 0x78, 0x45, 0x84, 0x10, 0x08, 0x09, 0xa6,
          0x09, 0x52, 0x6e, 0xaf, 0xd1, 0xa0, 0xf5, 0x4e, 0x61, 0x42, 0xe0, 0x29,
          0xc5, 0xc2, 0x1f, 0x3e, 0xf4, 0x15, 0x83, 0xbb, 0x77, 0xa8, 0x05, 0xef,
          0xbe, 0x1c, 0xe9, 0xd8, 0x62, 0x79, 0xe3, 0x26, 0x19, 0x67, 0x8e, 0xaa,
          0xb7, 0x85, 0x3f, 0xa0, 0xf0, 0x2a, 0x16, 0xab, 0xb6, 0x41, 0x66, 0xad,
          0x0c, 0xda, 0xc7, 0x53, 0x07, 0xeb, 0x38, 0x74, 0xb8, 0x4c, 0xb0, 0xa1,
          0x4a, 0xa5, 0xa2, 0x41, 0x08, 0xcc, 0x06, 0x88, 0xc4, 0x8c, 0x9a, 0x4a,
          0x1e, 0x82, 0x6d, 0xdb, 0xe6, 0xf0, 0xb1, 0x22, 0x19, 0x67, 0x8e, 0x8a,
          0xce, 0x61, 0x05, 0x02, 0xcc, 0xdf, 0xaf, 0xe3, 0xf6, 0x8d, 0x3c, 0x33,
          0xbf, 0x6c, 0x01, 0xc5, 0x9d, 0xcc, 0x20, 0x9d, 0x47, 0x2c, 0x4e, 0x0d,
          0x35, 0xd2, 0xf3, 0xba, 0x81, 0xd2, 0x05, 0xd0, 0x60, 0x05, 0xe4, 0xde,
          0x60, 0xc7, 0x71, 0x70, 0xbd, 0x22, 0x55, 0xaf, 0x80, 0x15, 0x08, 0x70,
          0xf7, 0x0e, 0x8c, 0x0e, 0xaf, 0x00, 0x92, 0x60, 0xa3, 0x49, 0x57, 0x77,
          0x14, 0xe9, 0x13, 0xa4, 0x17, 0x5d, 0x16, 0x92, 0x5b, 0x5c, 0xbe, 0x50,
          0x60, 0x79, 0x31, 0xcc, 0x3b, 0x67, 0x7c, 0x68, 0xca, 0x28, 0xa5, 0xf7,
          0x06, 0x6f, 0x5b, 0x26, 0xf0, 0x07, 0x60, 0xfe, 0xbe, 0xc9, 0xe8, 0xf0,
          0x32, 0xa0, 0xe8, 0x1b, 0x68, 0xa1, 0x7f, 0xc0, 0x20, 0x14, 0x15, 0x98,
          0x96, 0x47, 0xa9, 0xd8, 0xc4, 0xc4, 0x98, 0xcb, 0xe8, 0xf0, 0x5f, 0x5c,
          0xff, 0x3a, 0x4d, 0x4b, 0xac, 0x8b, 0xfe, 0x21, 0x0b, 0x29, 0x0b, 0x35,
          0xe0, 0x1a, 0xfd, 0x42, 0x82, 0x57, 0xb1, 0xb8, 0x75, 0xcd, 0x01, 0xaa,
          0xf4, 0xbd, 0x1d, 0xe5, 0xdd, 0x0f, 0xea, 0x68, 0x8e, 0x15, 0x30, 0xad,
          0x1c, 0x1b, 0xab, 0x16, 0xb7, 0xaf, 0x17, 0x39, 0x35, 0x58, 0xcf, 0x99,
          0x73, 0x11, 0xc0, 0xe4, 0x9b, 0xaf, 0xd6, 0x70, 0xd6, 0x25, 0x86, 0x21,
          0xf6, 0x57, 0x6c, 0x9a, 0xb0, 0x9a, 0x32, 0x48, 0xfe, 0x5a, 0xe6, 0xc0,
          0x81, 0x03, 0xf4, 0x9f, 0x36, 0x51, 0x72, 0x13, 0xd3, 0x30, 0x70, 0x36,
          0x9b, 0xf8, 0xe2, 0x52, 0x9e, 0xe9, 0xc9, 0x0c, 0xeb, 0x2b, 0x2e, 0xef,
          0x7f, 0x14, 0x61, 0x7c, 0x2c, 0x48, 0x72, 0x2a, 0x4b, 0x72, 0x26, 0x84,
          0x51, 0x43, 0xfa, 0x87, 0x62, 0x29, 0x05, 0xeb, 0xeb, 0x65, 0x14, 0x25,
          0x3a, 0x9e, 0xb5, 0x08, 0x45, 0x04, 0xa6, 0x69, 0x50, 0xc8, 0x86, 0x18,
          0xbe, 0xe0, 0x30, 0x3d, 0x99, 0x21, 0x7c, 0xb0, 0x9e, 0x9e, 0xd7, 0x1a,
          0x09, 0x36, 0x95, 0x38, 0xfa, 0x52, 0x3d, 0x00, 0xe9, 0xa5, 0x0a, 0x42,
          0xee, 0xd3, 0xbc, 0x6d, 0x8f, 0x41, 0x28, 0x1f, 0xa0, 0x90, 0x3e, 0x85,
          0x61, 0x79, 0x6c, 0xac, 0x86, 0xf8, 0xf2, 0xb3, 0x3c, 0x33, 0x93, 0x0e,
          0xa0, 0xf0, 0xfb, 0x4d, 0x5a, 0xdb, 0xeb, 0x30, 0xad, 0x22, 0x86, 0x6f,
          0xfb, 0xba, 0xd2, 0x26, 0xd2, 0xa8, 0x05, 0xd7, 0xee, 0x04, 0x84, 0x22,
          0x02, 0x08, 0x92, 0x5a, 0xaa, 0x50, 0xde, 0xaa, 0x67, 0xe2, 0x7b, 0x97,
          0xe9, 0xc9, 0x0c, 0xcd, 0xcd, 0xf5, 0xc4, 0x62, 0xcd, 0xd8, 0xa9, 0x4d,
          0x2e, 0x7e, 0x62, 0xf3, 0x60, 0x32, 0x44, 0xea, 0x41, 0x09, 0x89, 0xa4,
          0xed, 0x19, 0x1f, 0x52, 0xa8, 0xfd, 0x15, 0x6b, 0x05, 0x6d, 0x9d, 0x55,
          0x3a, 0x9e, 0x0b, 0xf0, 0xe7, 0xef, 0x45, 0xc6, 0xbf, 0x73, 0x39, 0xfd,
          0x5e, 0x90, 0xec, 0x7a, 0x88, 0x57, 0x5e, 0x6d, 0x20, 0xd6, 0xe6, 0xe7,
          0xe2, 0xa7, 0x55, 0xec, 0x85, 0x4d, 0xce, 0x7f, 0xb8, 0x88, 0x5b, 0x2c,
          0xd0, 0x14, 0xf6, 0xf3, 0xc2, 0x71, 0x4d, 0x3a, 0x55, 0x0b, 0xae, 0x55,
          0xac, 0x35, 0xc1, 0x06, 0x8f, 0xb7, 0x06, 0x1b, 0x00, 0xb8, 0x7a, 0x79,
          0x95, 0x89, 0x9b, 0x25, 0xce, 0x9e, 0x0b, 0x73, 0xe2, 0x4d, 0x8f, 0xae,
          0x97, 0xb3, 0x7c, 0x7c, 0xbe, 0x93, 0x70, 0xb4, 0x99, 0xfc, 0xe6, 0x06,
          0x6e, 0xb5, 0x4a, 0xdf, 0x60, 0x98, 0xf6, 0xee, 0x2a, 0x9a, 0x7f, 0xf9,
          0xc7, 0x08, 0x41, 0xa5, 0x02, 0x3d, 0x27, 0x0d, 0xec, 0xa5, 0x10, 0xd7,
          0x46, 0x53, 0x7c, 0x7e, 0x29, 0xcd, 0xf8, 0xad, 0x20, 0x47, 0x5f, 0xb4,
          0x30, 0x7c, 0x82, 0xd4, 0xfc, 0x32, 0x4e, 0xc6, 0x85, 0x9d, 0x66, 0xfd,
          0xfc, 0x63, 0x81, 0xe4, 0x40, 0x14, 0xb3, 0x6e, 0x9f, 0xe6, 0xe5, 0x72,
          0x39, 0x3c, 0xa5, 0x40, 0x68, 0xaa, 0x22, 0xcf, 0xd0, 0x59, 0x8b, 0x68,
          0x6b, 0x37, 0x37, 0xae, 0xae, 0x90, 0xbc, 0x97, 0x25, 0x79, 0x6f, 0x03,
          0x50, 0x80, 0x24, 0x14, 0xb6, 0xe8, 0x1b, 0x7c, 0x9a, 0x9f, 0x7e, 0x70,
          0x70, 0x8a, 0x2e, 0x6e, 0x59, 0x91, 0xcf, 0xee, 0xe3, 0xf1, 0xec, 0xec,
          0x2c, 0x0b, 0x0b, 0x3e, 0x8c, 0x06, 0x81, 0x42, 0xa2, 0xa8, 0xd2, 0x37,
          0x28, 0x39, 0x71, 0x32, 0xc6, 0xdc, 0x6f, 0x2e, 0x2b, 0x8b, 0x0a, 0x8d,
          0xa4, 0xb5, 0xdd, 0xc7, 0xf3, 0xc7, 0x35, 0xed, 0x47, 0x14, 0x6f, 0x0c,
          0xb4, 0x50, 0x29, 0x83, 0x69, 0x3a, 0x4c, 0x25, 0xca, 0x7b, 0x83, 0x13,
          0x89, 0x04, 0xfa, 0x8a, 0x41, 0x24, 0x66, 0xe0, 0x0f, 0x48, 0xb4, 0x27,
          0x10, 0x86, 0xc2, 0x90, 0x12, 0x63, 0x67, 0x6c, 0x4a, 0x69, 0x90, 0xb6,
          0x15, 0xb6, 0xad, 0xe1, 0x5b, 0x85, 0x69, 0x4a, 0xf2, 0x39, 0x98, 0x4a,
          0x94, 0x58, 0x5b, 0xf6, 0xf6, 0x06, 0xbb, 0xae, 0x5b, 0x33, 0x4f, 0x1f,
          0x37, 0xe4, 0xa3, 0x53, 0x1e, 0x03, 0x1c, 0x8f, 0xc7, 0xc5, 0xa3, 0x12,
          0xff, 0x17, 0xf8, 0x49, 0xc3, 0xe3, 0xf1, 0xb8, 0xf8, 0x1b, 0x8b, 0xe6,
          0x90, 0x0a, 0xca, 0x9b, 0x61, 0xc9, 0x00, 0x00, 0x00, 0x00, 0x49, 0x45,
          0x4e, 0x44, 0xae, 0x42, 0x60, 0x82, 0x00, 0x00, 0x00, 0x00, 0x00
        };

      switch (id)
        {
        case image0_ID: {QImage img; img.loadFromData(image0_data, sizeof(image0_data), "PNG"); return QPixmap::fromImage(img);}
        default: return QPixmap();
        } // switch
    } // icon
  };

namespace Ui
  {
  class CQProgressItemText: public Ui_CQProgressItemText {};
} // namespace Ui

class CQProgressItemText : public CQProgressItem, public Ui::CQProgressItemText
  {
    Q_OBJECT

  public:
    CQProgressItemText(QWidget* parent = 0, const char* name = 0);
    ~CQProgressItemText();

    virtual bool initFromProcessReportItem(CProcessReportItem * pItem);
    virtual bool process();
    virtual bool reset();

  protected slots:
    virtual void languageChange();

  private:
    void (CQProgressItemText::*mpSetValue)();
    CCopasiParameter::Value mParameterValue;

    virtual void setValueFromDOUBLE();
    virtual void setValueFromINT();
    virtual void setValueFromUINT();
    virtual void setValueFromSTRING();
  };

#endif // CQPROGRESSITEMTEXT_H
