﻿using CodeHelper;
using Controls;
using DataBaseLib;
using HardWares.温度控制器;
using HardWares.温度控制器.SRS_PTC10;
using HardWares.源表;
using HardWares.端口基类;
using HardWares.端口基类部分;
using HardWares.纳米位移台;
using HardWares.纳米位移台.PI;
using ODMR_Lab.Windows;
using ODMR_Lab.位移台部分;
using ODMR_Lab.基本控件;
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

namespace ODMR_Lab.源表部分
{
    /// <summary>
    /// Page1.xaml 的交互逻辑
    /// </summary>
    public partial class DevicePage : PageBase
    {

        public List<PowerMeterInfo> PowerMeterList { get; set; } = new List<PowerMeterInfo>();

        /// <summary>
        /// 当前选中的源表
        /// </summary>
        public PowerMeterInfo CurrentSelectedMeter { get; set; } = null;

        public DevicePage()
        {
            InitializeComponent();
        }

        public override void Init()
        {
            CreateMeasureThread();
        }

        public override void CloseBehaviour()
        {
            MeasureListener?.Abort();
        }

        public override void UpdateDataBaseToUI()
        {
            SamplePointBox.Text = PowerMeterSamplePoint.ToString();
            SampleTimeBox.Text = PowerMeterSampleTime.ToString();
        }

        /// <summary>
        /// 列出所有需要向数据库取回的数据
        /// </summary>
        public override void ListDataBaseData()
        {
            DataBase.RegistateVar("PowerMeterSampleTime", this);
            DataBase.RegistateVar("PowerMeterSamplePoint", this);
        }

        private void NewConnect(object sender, RoutedEventArgs e)
        {
            ConnectWindow window = new ConnectWindow(typeof(PowerSourceBase), MainWindow.Handle);
            bool res = window.ShowDialog(out PortObject dev);
            if (res == true)
            {
                PowerMeterInfo controller = (PowerMeterInfo)new PowerMeterInfo().CreateDeviceInfo(dev as PowerSourceBase, window.ConnectInfo);

                controller.MaxSavePoint = PowerMeterSamplePoint;

                PowerMeterList.Add(controller);

                RefreshPanels();
            }
            else
            {
                return;
            }
        }

        private DecoratedButton CreateMeterBar(PowerMeterInfo device)
        {
            DecoratedButton grid = new DecoratedButton();
            grid.Height = 40;
            Btn.CloneStyleTo(grid);
            grid.Text = device.Device.ProductName;
            grid.ToolTip = device.Device.ProductName;
            grid.Click += SelectBtn;
            grid.Tag = device;

            ContextMenu menu = new ContextMenu();
            menu.PanelMinWidth = 200;
            menu.ItemHeight = 40;
            DecoratedButton btn = new DecoratedButton()
            {
                Text = "关闭设备"
            };
            btn.Tag = device;
            Btn.CloneStyleTo(btn);
            btn.Click += CloseDevice;
            menu.Items.Add(btn);

            btn = new DecoratedButton()
            {
                Text = "参数设置"
            };
            btn.Tag = device;
            Btn.CloneStyleTo(btn);
            btn.Click += OpenParamWindow;
            menu.Items.Add(btn);

            menu.ApplyToControl(grid);

            return grid;
        }

        private void SelectBtn(object sender, RoutedEventArgs e)
        {
            DecoratedButton btn = sender as DecoratedButton;
            foreach (var item in DeviceListPanel.Children)
            {
                (item as DecoratedButton).KeepPressed = false;
            }
            btn.KeepPressed = true;
            CurrentSelectedMeter = btn.Tag as PowerMeterInfo;
        }

        public void RefreshPanels()
        {
            DeviceListPanel.Children.Clear();

            foreach (var item in PowerMeterList)
            {
                DecoratedButton g = CreateMeterBar(item);
                DeviceListPanel.Children.Add(g);
            }
        }

        private void OpenParamWindow(object sender, RoutedEventArgs e)
        {
            PowerMeterInfo obj = (sender as DecoratedButton).Tag as PowerMeterInfo;
            ParameterWindow window = new ParameterWindow(obj.Device.AvailableParameterNames());
            window.ShowDialog();
        }

        private void CloseDevice(object sender, RoutedEventArgs e)
        {
            if (MessageWindow.ShowMessageBox("提示", "确定要关闭此设备吗？", MessageBoxButton.YesNo, owner: MainWindow.Handle) == MessageBoxResult.Yes)
            {
                PowerMeterInfo mover = (sender as DecoratedButton).Tag as PowerMeterInfo;
                mover.CloseDeviceInfoAndSaveParams(out bool result);
                if (result == false)
                {
                    return;
                }
                PowerMeterList.Remove(mover);
                CurrentSelectedMeter = null;
                RefreshPanels();
            }
        }

        #region 电流监控线程

        /// <summary>
        /// 采样时间间隔
        /// </summary>
        public double PowerMeterSampleTime { get; set; } = 0.1;


        private int maxSamplePoint = 100000;
        /// <summary>
        /// 采样最大点数
        /// </summary>
        public int PowerMeterSamplePoint
        {
            get
            {
                return maxSamplePoint;
            }
            set
            {
                maxSamplePoint = value;
                foreach (var item in PowerMeterList)
                {
                    item.MaxSavePoint = value;
                }
            }
        }

        private void ApplySet(object sender, RoutedEventArgs e)
        {
            try
            {
                PowerMeterSamplePoint = int.Parse(SamplePointBox.Text);
                PowerMeterSampleTime = double.Parse(SampleTimeBox.Text);
            }
            catch (Exception ex)
            {
                MessageWindow.ShowTipWindow("参数格式错误", MainWindow.Handle);
            }
        }

        Thread MeasureListener = null;

        bool IsSampling = false;

        /// <summary>
        /// 创建监听线程
        /// </summary>
        private void CreateMeasureThread()
        {
            MeasureListener = new Thread(() =>
            {
                while (true)
                {
                    bool AllowSample = false;
                    Dispatcher.Invoke(() =>
                    {
                        AllowSample = AutoSample.IsSelected;
                    });
                    if (!AllowSample)
                    {
                        Dispatcher.Invoke(() =>
                        {
                            VoltagePanel.Value.Text = "--- V";
                            CurrentPanel.Value.Text = "--- A";
                        });
                        Thread.Sleep((int)(PowerMeterSampleTime * 1000));
                        continue;
                    }
                    IsSampling = true;
                    for (int i = 0; i < PowerMeterList.Count; i++)
                    {
                        if (PowerMeterList[i].AllowAutoMeasure)
                        {
                            PowerMeterList[i].IsMeasuring = true;
                            MeasureGroup result = PowerMeterList[i].Device.Measure();
                            PowerMeterList[i].IsMeasuring = false;
                            PowerMeterList[i].AddMeasurePoint(result.TimeStamp, result.Voltage, result.Current);
                            //更新显示面板
                            if (PowerMeterList[i] == CurrentSelectedMeter)
                            {
                                Dispatcher.Invoke(() =>
                                {
                                    VoltagePanel.Value.Text = result.Voltage.ToString() + " V";
                                    CurrentPanel.Value.Text = result.Current.ToString() + " A";
                                    CurrentPanel.SetRestrictState(result.IsCurrentLimited);
                                });
                            }
                        }
                    }
                    if (CurrentSelectedMeter == null)
                    {
                        Dispatcher.Invoke(() =>
                        {
                            VoltagePanel.Value.Text = "--- V";
                            CurrentPanel.Value.Text = "--- A";
                        });
                    }
                    IsSampling = false;
                    Thread.Sleep((int)(PowerMeterSampleTime * 1000));
                }
            });
            MeasureListener.Start();
        }

        #endregion


        ChartViewerWindow GraphWindow = new ChartViewerWindow(true);
        /// <summary>
        /// 显示图表
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void ShowGraph(object sender, RoutedEventArgs e)
        {
            if (CurrentSelectedMeter == null) return;
            List<ChartData1D> list = new List<ChartData1D>();
            list.Add(new NumricChartData1D() { Data = new List<double>(CurrentSelectedMeter.CurrentBuffer.ToArray()), Name = "电流值(A)" });
            list.Add(new NumricChartData1D() { Data = new List<double>(CurrentSelectedMeter.VoltageBuffer), Name = "电压值(V)" });
            list.Add(new TimeChartData1D() { Data = new List<DateTime>(CurrentSelectedMeter.Times), Name = "采样时间" });
            GraphWindow.Title = "源表数据：" + CurrentSelectedMeter.Device.ProductName;
            GraphWindow.ShowAs1D(list);
            GraphWindow.Topmost = true;
            GraphWindow.Topmost = false;
        }
    }
}
