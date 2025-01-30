﻿using CodeHelper;
using Controls;
using HardWares.温度控制器;
using HardWares.温度控制器.SRS_PTC10;
using HardWares.端口基类;
using HardWares.端口基类部分;
using HardWares.纳米位移台;
using HardWares.纳米位移台.PI;
using ODMR_Lab.Windows;
using ODMR_Lab.基本窗口;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using 温度监控程序.Windows;
using ContextMenu = Controls.ContextMenu;

namespace ODMR_Lab.位移台部分
{
    /// <summary>
    /// Page1.xaml 的交互逻辑
    /// </summary>
    public partial class DevicePage : PageBase
    {

        public List<NanoMoverInfo> MoverList { get; set; } = new List<NanoMoverInfo>();


        public DevicePage()
        {
            InitializeComponent();
        }

        public override void Init()
        {
        }

        public override void CloseBehaviour()
        {
        }

        private void NewConnect(object sender, RoutedEventArgs e)
        {
            ConnectWindow window = new ConnectWindow(typeof(NanoControllerBase), MainWindow.Handle);
            bool res = window.ShowDialog(out PortObject dev);
            if (res == true)
            {
                DecoratedButton btn = sender as DecoratedButton;

                NanoMoverInfo controller = new NanoMoverInfo() { Device = dev as NanoControllerBase, ConnectInfo = window.ConnectInfo };
                controller.CreateDeviceInfoBehaviour();
                if (btn.Name == "MagnetBtn")
                {
                    foreach (var item in controller.Stages)
                    {
                        item.PartType = PartTypes.Magnnet;
                    }
                }
                if (btn.Name == "SampleBtn")
                {
                    foreach (var item in controller.Stages)
                    {
                        item.PartType = PartTypes.Sample;
                    }
                }
                if (btn.Name == "MicrowaveBtn")
                {
                    foreach (var item in controller.Stages)
                    {
                        item.PartType = PartTypes.Microwave;
                    }
                }
                if (btn.Name == "ProbeBtn")
                {
                    foreach (var item in controller.Stages)
                    {
                        item.PartType = PartTypes.Probe;
                    }
                }
                MoverList.Add(controller);
                RefreshPanels();
            }
            else
            {
                return;
            }
            NanoControllerBase tem = dev as NanoControllerBase;
        }


        public void RefreshPanels()
        {
            ProbeMoverList.ClearItems();
            SampleMoverList.ClearItems();
            MWMoverList.ClearItems();
            MagnetMoverList.ClearItems();

            foreach (var item in MoverList)
            {
                foreach (var stage in item.Stages)
                {
                    if (stage.PartType == PartTypes.Probe)
                    {
                        ProbeMoverList.AddItem(stage, stage.Parent.Device.ProductName, stage.Device.AxisName);
                    }
                    if (stage.PartType == PartTypes.Sample)
                    {
                        SampleMoverList.AddItem(stage, stage.Parent.Device.ProductName, stage.Device.AxisName);
                    }
                    if (stage.PartType == PartTypes.Magnnet)
                    {
                        MagnetMoverList.AddItem(stage, stage.Parent.Device.ProductName, stage.Device.AxisName);
                    }
                    if (stage.PartType == PartTypes.Microwave)
                    {
                        MWMoverList.AddItem(stage, stage.Parent.Device.ProductName, stage.Device.AxisName);
                    }
                }
            }
        }


        private void OpenLabelWindow(object sender, RoutedEventArgs e)
        {
            DecoratedButton btn = sender as DecoratedButton;
            try
            {
                MoverTestWindow window = null;
                if (btn.Name == "MagnetLabelBtn")
                {
                    window = new MoverTestWindow(PartTypes.Magnnet, "磁铁位移台设置窗口");
                }
                if (btn.Name == "PrbeLabelBtn")
                {
                    window = new MoverTestWindow(PartTypes.Probe, "探针位移台设置窗口");
                }
                if (btn.Name == "SampleLabelBtn")
                {
                    window = new MoverTestWindow(PartTypes.Sample, "样品位移台设置窗口");
                }
                if (btn.Name == "MicrowaveLabelBtn")
                {
                    window = new MoverTestWindow(PartTypes.Microwave, "微波线位移台设置窗口");
                }

                window?.ShowDialog();
            }
            catch (Exception ex)
            {
                MessageWindow.ShowTipWindow("设备被占用", MainWindow.Handle);
                return;
            }
        }

        #region 设备列表
        private void ContextMenuEvent(int arg1, int arg2, object arg3)
        {
            NanoStageInfo info = arg3 as NanoStageInfo;
            #region 关闭设备
            if (arg1 == 0)
            {
                if (MessageWindow.ShowMessageBox("提示", "确定要关闭此设备吗？", MessageBoxButton.YesNo, owner: MainWindow.Handle) == MessageBoxResult.Yes)
                {
                    info.Parent.CloseDeviceInfoAndSaveParams(out bool result);
                    if (result == false) return;
                    MoverList.Remove(info.Parent);
                    RefreshPanels();
                }
            }
            #endregion
            #region 参数设置
            if (arg1 == 1)
            {
                ParameterWindow window = new ParameterWindow(info.Device.AvailableParameterNames());
                window.ShowDialog();
            }
            #endregion
        }
        #endregion
    }
}
